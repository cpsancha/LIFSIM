/*
 * mainV03_56.c
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

#include "mainV03_56.h"
#include "mainV03_56_private.h"

/* Block signals (auto storage) */
B_mainV03_56_T mainV03_56_B;

/* Continuous states */
X_mainV03_56_T mainV03_56_X;

/* Block states (auto storage) */
DW_mainV03_56_T mainV03_56_DW;

/* Previous zero-crossings (trigger) states */
PrevZCX_mainV03_56_T mainV03_56_PrevZCX;

/* Real-time model */
RT_MODEL_mainV03_56_T mainV03_56_M_;
RT_MODEL_mainV03_56_T *const mainV03_56_M = &mainV03_56_M_;
static void rate_scheduler(void);

/*     Initialize pressure and temperature tables. */
void InitCalcAtmosCOESA(real_T *temperature76, real_T *pressureRatio76)
{
  if (temperature76[0] != TEMPERATURE0 ) {
    int_T k;
    temperature76[0] = TEMPERATURE0;
    pressureRatio76[0] = 1.0;

    /* set up the data at the 1976 altitude breakpoints */
    for (k=0; k<(NUM1976PTS-1); k++) {
      if (tempGradient76[k] != 0.0) {
        temperature76[k+1] = temperature76[k] +
          tempGradient76[k]*(altitude76[k+1] - altitude76[k]);
        pressureRatio76[k+1] = pressureRatio76[k] *
          exp(log(temperature76[k]/temperature76[k+1]) * GMR/tempGradient76[k]);
      } else {
        temperature76[k+1] = temperature76[k];
        pressureRatio76[k+1] = pressureRatio76[k] *
          exp((-GMR)*(altitude76[k+1] - altitude76[k])/temperature76[k]);
      }
    }
  }
}

/*
 *     Using cached pressure and temperature tables, find the
 *     working interval and perform logarithmic interpolation.
 */
void CalcAtmosCOESA(const real_T *altitude, real_T *temp, real_T *pressure,
                    real_T *density, real_T *speedofsound, real_T *temperature76,
                    real_T *pressureRatio76, int_T numPoints)
{
  int_T i;
  for (i=0; i < numPoints; i++) {
    int_T bottom = 0;
    int_T top = NUM1976PTS-1;
    int_T idx;

    /* Find altitude interval using binary search
     *
     * Deal with the extreme cases first:
     *   if altitude <= altitude76[bottom] then return idx = bottom
     *   if altitude >= altitude76[top]    then return idx = top
     */
    if (altitude[i] <= altitude76[bottom]) {
      idx = bottom;
    } else if (altitude[i] >= altitude76[top]) {
      idx = NUM1976PTS-2;
    } else {
      for (;;) {
        idx = (bottom + top)/2;
        if (altitude[i] < altitude76[idx]) {
          top = idx - 1;
        } else if (altitude[i] >= altitude76[idx+1]) {
          bottom = idx + 1;
        } else {
          /* we have altitude76[idx] <= altitude[i] < altitude76[idx+1],
           * so break and just use idx
           */
          break;
        }
      }
    }

    /* Interval has been obtained, now do linear temperature
     * interpolation and log pressure interpolation.
     */
    if (tempGradient76[idx] != 0.0 ) {
      temp[i] = temperature76[idx] +
        tempGradient76[idx] * (altitude[i] - altitude76[idx]);
      pressure[i] = PRESSURE0 * pressureRatio76[idx] *
        (rt_powd_snf(temperature76[idx]/temp[i], GMR/tempGradient76[idx]));
    } else {
      temp[i] = temperature76[idx];
      pressure[i] = PRESSURE0 * pressureRatio76[idx] *
        exp((-GMR)*(altitude[i] - altitude76[idx]) / temperature76[idx]);
    }

    density[i] = pressure[i] / ((R_HAT/MOL_WT)*temp[i]);
    speedofsound[i] = sqrt(GAMMA*temp[i]*(R_HAT/MOL_WT));
  }
}

/*     Calculate the Julian date */
real_T calc_Julian_date(real_T uDataMonth, real_T uDataDay, real_T uDataYear)
{
  real_T month, day, year;
  real_T temp1, temp2;
  month = uDataMonth;
  day = uDataDay;
  year = uDataYear;

  /* calculate Julian Date (JD) */
  if ((MonthIdx)uDataMonth <= FEBRUARY ) {
    year -= 1.0;
    month += 12.0;
  }

  temp1 = floor(year/100.0);
  temp2 = 2.0 - temp1 + floor(temp1/4.0);
  return floor(365.25*(year + 4716.0)) +
    floor(30.6001*( month + 1.0)) + temp2 + day - 1524.0;
}

/*     Calculate the variables shared by close approx. and exact methods. */
int_T wgs84_calc_shared_vars(real_T uDataMonth, real_T uDataDay, real_T
  uDataYear, int_T uDataPrecessing, int_T uDataCentrifugal, real_T h, real_T E2,
  real_T cosphi, real_T sinphi, real_T sin2phi, real_T coslambda, real_T
  sinlambda, real_T GM, real_T *gamma_u_ptr, real_T *gamma_beta_ptr, real_T
  *cosbeta_ptr, real_T *sinbeta_ptr, real_T *u_ptr, real_T *u2E2_ptr, real_T
  *w_ptr)
{
  real_T N, x_rec, y_rec, z_rec, D, u2, omega;
  real_T beta, sin2beta, cos2beta, q, qo, q_prime, cf_u, cf_beta;
  real_T w, u, u2E2, sinbeta, cosbeta, gamma_u, gamma_beta;

  /* Radius of Curvature in prime vertical (N) /eq. 4-15/ */
  N = (WGS84_A)/(sqrt(1.0 - (WGS84_E_2)*sin2phi));

  /* Calculate rectangular coordinates /eq. 4-14/ */
  x_rec = ( N + h )*cosphi*coslambda;
  y_rec = ( N + h )*cosphi*sinlambda;
  z_rec = ((WGS84_B_A)*(WGS84_B_A)*N + h )*sinphi;

  /* Calculate various parameters */
  D = x_rec*x_rec + y_rec*y_rec + z_rec*z_rec - E2;
  u2 = 0.5*D*( 1.0 + sqrt(1.0 + 4.0*E2*z_rec*z_rec/(D*D)));
  u2E2 = u2 + E2;

  /* /eq. 4-8/ */
  u = sqrt(u2);

  /* /eq. 4-9/ */
  beta = atan(z_rec*sqrt(u2E2)/(u*sqrt(x_rec*x_rec + y_rec*y_rec)));

  /* generate common sines and cosines */
  sinbeta = sin(beta);
  sin2beta = sinbeta*sinbeta;
  cosbeta = cos(beta);
  cos2beta = cosbeta*cosbeta;

  /* /eq. 4-10/ */
  w = sqrt(( u2 + E2*sin2beta )/( u2E2 ));

  /* /eq. 4-11/ */
  q = 0.5*(( 1.0 + 3.0*u2/( E2 ))*atan((WGS84_EL)/u) - 3.0*u/(WGS84_EL));

  /* /eq. 4-12/ */
  qo = 0.5*(( 1.0 + 3.0*(WGS84_B)*(WGS84_B)/( E2 ))*atan((WGS84_EL)/(WGS84_B)) -
            3.0*(WGS84_B)/(WGS84_EL));

  /* /eq. 4-13/ */
  q_prime = 3.0*(( 1.0 + u2/( E2 ))*( 1.0 - (u/(WGS84_EL))*atan((WGS84_EL)/u)))
    - 1.0;

  /* Use precessing reference frame? */
  if (uDataPrecessing == 0 ) {
    omega = WGS84_W_DEF;
  } else {
    real_T JD;
    if (uDataYear == 99.0 ) {
      JD = uDataDay;
    } else {
      /* calculate Julian Date (JD) */
      JD = calc_Julian_date( uDataMonth, uDataDay, uDataYear );
    }

    /* Ang Vel of Earth (rad/sec) [precessing ref. frame] /eq. 3-8/ */
    omega = WGS84_W_PRM + 7.086e-12 +
      4.3e-15*((JD - 2451545.0)/36525.0);
  }

  /* Use Centrifugal Force? */
  if (uDataCentrifugal == 0 ) {
    cf_u = u*cos2beta*omega*omega/w;
    cf_beta = sqrt(u2E2)*cosbeta*sinbeta*omega*omega/w;
  } else {
    cf_u = 0.0;
    cf_beta = 0.0;
  }

  /* /eq. 4-5/ */
  gamma_u = -( GM/u2E2 + omega*omega*(WGS84_A)*(WGS84_A)*(WGS84_EL)*q_prime*
              ( 0.5*sin2beta - 1.0/6.0 )/(u2E2*qo))/w + cf_u;

  /* /eq. 4-6/ */
  gamma_beta = omega*omega*(WGS84_A)*(WGS84_A)*q*sinbeta*cosbeta/
    (sqrt(u2E2)*w*qo) - cf_beta;
  *gamma_u_ptr = gamma_u;
  *gamma_beta_ptr = gamma_beta;
  *cosbeta_ptr = cosbeta;
  *sinbeta_ptr = sinbeta;
  *u_ptr = u;
  *u2E2_ptr = u2E2;
  *w_ptr = w;
  return 0;
}

/*     Exact WGS84 model of ellipsoid normal gravity */
void wgs84_exact(real_T *h, real_T *phi, real_T *lambda, real_T uDataMonth,
                 real_T uDataDay, real_T uDataYear, int_T uDataPrecessing, int_T
                 uDataCentrifugal, real_T E2, real_T GM, real_T opt_m2ft, real_T
                 *y, real_T *gamma_h, real_T *gamma_phi,int_T k)
{
  /* Calculating exact normal gravity */
  real_T psi, tanphi, sinpsi, cospsi, alpha, cosalpha, sinalpha;
  real_T gamma_r, gamma_psi;
  real_T gamma_u, gamma_beta, u, u2E2, sinbeta, cosbeta, w;
  real_T sinphi, sin2phi, cosphi, coslambda, sinlambda;
  int_T i;
  for (i = 0; i < k; i++ ) {
    /* generate common sines and cosines of lat and long angles */
    sinphi = sin(phi[i]);
    sin2phi = sinphi*sinphi;
    cosphi = cos(phi[i]);
    coslambda = cos(lambda[i]);
    sinlambda = sin(lambda[i]);
    tanphi = sinphi/cosphi;

    /* Calc. geocentric latitude (psi) from geodetic latitude (phi) */
    psi = atan(tanphi*( 1.0 - 1.0/(WGS84_INV_F))*( 1.0 - 1.0/(WGS84_INV_F)));
    cospsi = cos(psi);
    sinpsi = sin(psi);

    /* /eq. 4-20/ */
    alpha = phi[i] - psi;
    cosalpha = cos(alpha);
    sinalpha = sin(alpha);
    wgs84_calc_shared_vars( uDataMonth,uDataDay,uDataYear,uDataPrecessing,
      uDataCentrifugal,h[i],E2,cosphi,sinphi,sin2phi,coslambda,
      sinlambda,GM,&gamma_u,&gamma_beta,&cosbeta,
      &sinbeta,&u,&u2E2,&w );

    /* /eq. 4-17, 4-18, 4-19/
     *
     * |gamma_r     |         |gamma_u     |
     * |gamma_psi   | = R2*R1*|gamma_beta  |
     * |gamma_lambda|         |gamma_lambda|
     *
     * where:
     * gamma_lambda = 0
     *
     *      |cosbeta*coslambda*u   -sinbeta*coslambda    -sinlambda|
     *      |-------------------   ------------------              |
     *      |    w*sqrt(u2E2)              w                       |
     * R1 = |                                                      |
     *      |cosbeta*sinlambda*u   -sinbeta*sinlambda     coslambda|
     *      |-------------------   ------------------              |
     *      |    w*sqrt(u2E2)              w                       |
     *      |                                                      |
     *      |      sinbeta              cosbeta*u            0     |
     *      |-------------------   ------------------              |
     *      |        w                w*sqrt(u2E2)                 |
     *
     *
     *      | cospsi*coslambda  cospsi*sinlambda  sinpsi|
     * R2 = |-sinpsi*coslambda -sinpsi*sinlambda  cospsi|
     *      |-sinlambda         coslambda           0   |
     *
     *
     *         |cosB*cosP*u    sinB*sinP  sinP*cosB*u    sinB*cosP   0|
     *         |----------- +  ---------  ----------- -  ---------    |
     *         |w*sqrt(u2E2)       w      w*sqrt(u2E2)       w        |
     * R2*R1 = |                                                      |
     *         |sinB*cosP   sinP*cosB*u  cosB*cosP*u   sinB*sinP     0|
     *         |--------- - -----------  ----------- + ----------     |
     *         |    w       w*sqrt(u2E2) w*sqrt(u2E2)      w          |
     *         |                                                      |
     *         |           0                         0               1|
     *
     *
     *        where cosB = cosbeta
     *              sinB = sinbeta
     *              cosP = cospsi
     *              sinP = sinpsi
     */
    gamma_r = (cosbeta*cospsi*u/(w*sqrt(u2E2))+sinbeta*sinpsi/w)*
      gamma_u + (sinpsi*cosbeta*u/(w*sqrt(u2E2))-
                 sinbeta*cospsi/w)*gamma_beta;
    gamma_psi = (sinbeta*cospsi/w-sinpsi*cosbeta*u/(w*sqrt(u2E2)))*
      gamma_u + (cosbeta*cospsi*u/(w*sqrt(u2E2))+
                 sinbeta*sinpsi/w)*gamma_beta;

    /* "normal" gravity  / eq. 4-16 / (positive downwards)*/
    gamma_h[i] = (-gamma_r)*cosalpha - gamma_psi*sinalpha;

    /* "tangent" gravity  / eq. 4-23 / (positive northward)*/
    gamma_phi[i] = (-gamma_r)*sinalpha + gamma_psi*cosalpha;

    /* Return total normal gravity as the output / eq. 4-24 / */
    y[i] = opt_m2ft*sqrt(gamma_h[i]*gamma_h[i] + gamma_phi[i]*gamma_phi[i]);
  }
}

int32_T plook_s32dd_bincp(real_T u, const real_T bp[], uint32_T maxIndex, real_T
  *fraction, int32_T *prevIndex)
{
  int32_T bpIndex;

  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Clip'
     Use previous index: 'on'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u <= bp[0U]) {
    bpIndex = 0;
    *fraction = 0.0;
  } else if (u < bp[maxIndex]) {
    bpIndex = binsearch_s32d_prevIdx(u, bp, (uint32_T)*prevIndex, maxIndex);
    *fraction = (u - bp[(uint32_T)bpIndex]) / (bp[bpIndex + 1U] - bp[(uint32_T)
      bpIndex]);
  } else {
    bpIndex = (int32_T)(maxIndex - 1U);
    *fraction = 1.0;
  }

  *prevIndex = bpIndex;
  return bpIndex;
}

uint32_T plook_bincpa(real_T u, const real_T bp[], uint32_T maxIndex, real_T
                      *fraction, uint32_T *prevIndex)
{
  uint32_T bpIndex;

  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Clip'
     Use previous index: 'on'
     Use last breakpoint for index at or above upper limit: 'on'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u <= bp[0U]) {
    bpIndex = 0U;
    *fraction = 0.0;
  } else if (u < bp[maxIndex]) {
    bpIndex = binsearch_u32d_prevIdx(u, bp, *prevIndex, maxIndex);
    *fraction = (u - bp[bpIndex]) / (bp[bpIndex + 1U] - bp[bpIndex]);
  } else {
    bpIndex = maxIndex;
    *fraction = 0.0;
  }

  *prevIndex = bpIndex;
  return bpIndex;
}

real_T intrp2d_la_pw(const uint32_T bpIndex[], const real_T frac[], const real_T
                     table[], uint32_T stride, const uint32_T maxIndex[])
{
  real_T y;
  real_T yR_1d;
  uint32_T offset_1d;

  /* Interpolation 2-D
     Interpolation method: 'Linear'
     Use last breakpoint for index at or above upper limit: 'on'
     Overflow mode: 'portable wrapping'
   */
  offset_1d = bpIndex[1U] * stride + bpIndex[0U];
  if (bpIndex[0U] == maxIndex[0U]) {
    y = table[offset_1d];
  } else {
    y = (table[offset_1d + 1U] - table[offset_1d]) * frac[0U] + table[offset_1d];
  }

  if (bpIndex[1U] == maxIndex[1U]) {
  } else {
    offset_1d += stride;
    if (bpIndex[0U] == maxIndex[0U]) {
      yR_1d = table[offset_1d];
    } else {
      yR_1d = (table[offset_1d + 1U] - table[offset_1d]) * frac[0U] +
        table[offset_1d];
    }

    y += (yR_1d - y) * frac[1U];
  }

  return y;
}

int32_T plook_s32dd_bincpa(real_T u, const real_T bp[], uint32_T maxIndex,
  real_T *fraction, int32_T *prevIndex)
{
  int32_T bpIndex;

  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Clip'
     Use previous index: 'on'
     Use last breakpoint for index at or above upper limit: 'on'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u <= bp[0U]) {
    bpIndex = 0;
    *fraction = 0.0;
  } else if (u < bp[maxIndex]) {
    bpIndex = binsearch_s32d_prevIdx(u, bp, (uint32_T)*prevIndex, maxIndex);
    *fraction = (u - bp[(uint32_T)bpIndex]) / (bp[bpIndex + 1U] - bp[(uint32_T)
      bpIndex]);
  } else {
    bpIndex = (int32_T)maxIndex;
    *fraction = 0.0;
  }

  *prevIndex = bpIndex;
  return bpIndex;
}

real_T intrp8d_s32dla_pw(const int32_T bpIndex[], const real_T frac[], const
  real_T table[], const uint32_T stride[], const uint32_T maxIndex[])
{
  real_T y;
  real_T yR;
  uint32_T offset;

  /* Interpolation 8-D
     Interpolation method: 'Linear'
     Use last breakpoint for index at or above upper limit: 'on'
     Overflow mode: 'portable wrapping'
   */
  offset = (uint32_T)(bpIndex[7U] * (int32_T)stride[7U]);
  y = intrp7d_s32dla_pw(bpIndex, frac, &table[offset], stride, maxIndex);
  if ((uint32_T)bpIndex[7U] == maxIndex[7U]) {
  } else {
    offset += stride[7U];
    yR = intrp7d_s32dla_pw(bpIndex, frac, &table[offset], stride, maxIndex);
    y += (yR - y) * frac[7U];
  }

  return y;
}

real_T intrp7d_s32dla_pw(const int32_T bpIndex[], const real_T frac[], const
  real_T table[], const uint32_T stride[], const uint32_T maxIndex[])
{
  real_T y;
  real_T yR;
  uint32_T offset;

  /* Interpolation 7-D
     Interpolation method: 'Linear'
     Use last breakpoint for index at or above upper limit: 'on'
     Overflow mode: 'portable wrapping'
   */
  offset = (uint32_T)(bpIndex[6U] * (int32_T)stride[6U]);
  y = intrp6d_s32dla_pw(bpIndex, frac, &table[offset], stride, maxIndex);
  if ((uint32_T)bpIndex[6U] == maxIndex[6U]) {
  } else {
    offset += stride[6U];
    yR = intrp6d_s32dla_pw(bpIndex, frac, &table[offset], stride, maxIndex);
    y += (yR - y) * frac[6U];
  }

  return y;
}

real_T intrp6d_s32dla_pw(const int32_T bpIndex[], const real_T frac[], const
  real_T table[], const uint32_T stride[], const uint32_T maxIndex[])
{
  real_T y;
  real_T yR;
  uint32_T offset;

  /* Interpolation 6-D
     Interpolation method: 'Linear'
     Use last breakpoint for index at or above upper limit: 'on'
     Overflow mode: 'portable wrapping'
   */
  offset = (uint32_T)(bpIndex[5U] * (int32_T)stride[5U]);
  y = intrp5d_s32dla_pw(bpIndex, frac, &table[offset], stride, maxIndex);
  if ((uint32_T)bpIndex[5U] == maxIndex[5U]) {
  } else {
    offset += stride[5U];
    yR = intrp5d_s32dla_pw(bpIndex, frac, &table[offset], stride, maxIndex);
    y += (yR - y) * frac[5U];
  }

  return y;
}

real_T intrp5d_s32dla_pw(const int32_T bpIndex[], const real_T frac[], const
  real_T table[], const uint32_T stride[], const uint32_T maxIndex[])
{
  real_T y;
  uint32_T offset_4d;
  real_T yR_1d;
  uint32_T offset_0d;
  real_T yL_1d;
  uint32_T offset_1d;
  real_T yL_1d_0;
  real_T yL_1d_1;

  /* Interpolation 5-D
     Interpolation method: 'Linear'
     Use last breakpoint for index at or above upper limit: 'on'
     Overflow mode: 'portable wrapping'
   */
  offset_4d = ((((uint32_T)(bpIndex[4U] * (int32_T)stride[4U]) + bpIndex[3U] *
                 (int32_T)stride[3U]) + bpIndex[2U] * (int32_T)stride[2U]) +
               bpIndex[1U] * (int32_T)stride[1U]) + bpIndex[0U];
  if ((uint32_T)bpIndex[0U] == maxIndex[0U]) {
    y = table[offset_4d];
  } else {
    y = (table[offset_4d + 1U] - table[offset_4d]) * frac[0U] + table[offset_4d];
  }

  if ((uint32_T)bpIndex[1U] == maxIndex[1U]) {
  } else {
    offset_0d = offset_4d + stride[1U];
    if ((uint32_T)bpIndex[0U] == maxIndex[0U]) {
      yR_1d = table[offset_0d];
    } else {
      yR_1d = (table[offset_0d + 1U] - table[offset_0d]) * frac[0U] +
        table[offset_0d];
    }

    y += (yR_1d - y) * frac[1U];
  }

  if ((uint32_T)bpIndex[2U] == maxIndex[2U]) {
  } else {
    offset_1d = offset_4d + stride[2U];
    if ((uint32_T)bpIndex[0U] == maxIndex[0U]) {
      yL_1d = table[offset_1d];
    } else {
      yL_1d = (table[offset_1d + 1U] - table[offset_1d]) * frac[0U] +
        table[offset_1d];
    }

    if ((uint32_T)bpIndex[1U] == maxIndex[1U]) {
    } else {
      offset_0d = offset_1d + stride[1U];
      if ((uint32_T)bpIndex[0U] == maxIndex[0U]) {
        yR_1d = table[offset_0d];
      } else {
        yR_1d = (table[offset_0d + 1U] - table[offset_0d]) * frac[0U] +
          table[offset_0d];
      }

      yL_1d += (yR_1d - yL_1d) * frac[1U];
    }

    y += (yL_1d - y) * frac[2U];
  }

  if ((uint32_T)bpIndex[3U] == maxIndex[3U]) {
  } else {
    offset_1d = offset_4d + stride[3U];
    if ((uint32_T)bpIndex[0U] == maxIndex[0U]) {
      yL_1d = table[offset_1d];
    } else {
      yL_1d = (table[offset_1d + 1U] - table[offset_1d]) * frac[0U] +
        table[offset_1d];
    }

    if ((uint32_T)bpIndex[1U] == maxIndex[1U]) {
    } else {
      offset_0d = offset_1d + stride[1U];
      if ((uint32_T)bpIndex[0U] == maxIndex[0U]) {
        yR_1d = table[offset_0d];
      } else {
        yR_1d = (table[offset_0d + 1U] - table[offset_0d]) * frac[0U] +
          table[offset_0d];
      }

      yL_1d += (yR_1d - yL_1d) * frac[1U];
    }

    if ((uint32_T)bpIndex[2U] == maxIndex[2U]) {
    } else {
      offset_1d += stride[2U];
      if ((uint32_T)bpIndex[0U] == maxIndex[0U]) {
        yL_1d_0 = table[offset_1d];
      } else {
        yL_1d_0 = (table[offset_1d + 1U] - table[offset_1d]) * frac[0U] +
          table[offset_1d];
      }

      if ((uint32_T)bpIndex[1U] == maxIndex[1U]) {
      } else {
        offset_0d = offset_1d + stride[1U];
        if ((uint32_T)bpIndex[0U] == maxIndex[0U]) {
          yR_1d = table[offset_0d];
        } else {
          yR_1d = (table[offset_0d + 1U] - table[offset_0d]) * frac[0U] +
            table[offset_0d];
        }

        yL_1d_0 += (yR_1d - yL_1d_0) * frac[1U];
      }

      yL_1d += (yL_1d_0 - yL_1d) * frac[2U];
    }

    y += (yL_1d - y) * frac[3U];
  }

  if ((uint32_T)bpIndex[4U] == maxIndex[4U]) {
  } else {
    offset_4d += stride[4U];
    if ((uint32_T)bpIndex[0U] == maxIndex[0U]) {
      yL_1d = table[offset_4d];
    } else {
      yL_1d = (table[offset_4d + 1U] - table[offset_4d]) * frac[0U] +
        table[offset_4d];
    }

    if ((uint32_T)bpIndex[1U] == maxIndex[1U]) {
    } else {
      offset_0d = offset_4d + stride[1U];
      if ((uint32_T)bpIndex[0U] == maxIndex[0U]) {
        yR_1d = table[offset_0d];
      } else {
        yR_1d = (table[offset_0d + 1U] - table[offset_0d]) * frac[0U] +
          table[offset_0d];
      }

      yL_1d += (yR_1d - yL_1d) * frac[1U];
    }

    if ((uint32_T)bpIndex[2U] == maxIndex[2U]) {
    } else {
      offset_1d = offset_4d + stride[2U];
      if ((uint32_T)bpIndex[0U] == maxIndex[0U]) {
        yL_1d_0 = table[offset_1d];
      } else {
        yL_1d_0 = (table[offset_1d + 1U] - table[offset_1d]) * frac[0U] +
          table[offset_1d];
      }

      if ((uint32_T)bpIndex[1U] == maxIndex[1U]) {
      } else {
        offset_0d = offset_1d + stride[1U];
        if ((uint32_T)bpIndex[0U] == maxIndex[0U]) {
          yR_1d = table[offset_0d];
        } else {
          yR_1d = (table[offset_0d + 1U] - table[offset_0d]) * frac[0U] +
            table[offset_0d];
        }

        yL_1d_0 += (yR_1d - yL_1d_0) * frac[1U];
      }

      yL_1d += (yL_1d_0 - yL_1d) * frac[2U];
    }

    if ((uint32_T)bpIndex[3U] == maxIndex[3U]) {
    } else {
      offset_1d = offset_4d + stride[3U];
      if ((uint32_T)bpIndex[0U] == maxIndex[0U]) {
        yL_1d_0 = table[offset_1d];
      } else {
        yL_1d_0 = (table[offset_1d + 1U] - table[offset_1d]) * frac[0U] +
          table[offset_1d];
      }

      if ((uint32_T)bpIndex[1U] == maxIndex[1U]) {
      } else {
        offset_0d = offset_1d + stride[1U];
        if ((uint32_T)bpIndex[0U] == maxIndex[0U]) {
          yR_1d = table[offset_0d];
        } else {
          yR_1d = (table[offset_0d + 1U] - table[offset_0d]) * frac[0U] +
            table[offset_0d];
        }

        yL_1d_0 += (yR_1d - yL_1d_0) * frac[1U];
      }

      if ((uint32_T)bpIndex[2U] == maxIndex[2U]) {
      } else {
        offset_1d += stride[2U];
        if ((uint32_T)bpIndex[0U] == maxIndex[0U]) {
          yL_1d_1 = table[offset_1d];
        } else {
          yL_1d_1 = (table[offset_1d + 1U] - table[offset_1d]) * frac[0U] +
            table[offset_1d];
        }

        if ((uint32_T)bpIndex[1U] == maxIndex[1U]) {
        } else {
          offset_0d = offset_1d + stride[1U];
          if ((uint32_T)bpIndex[0U] == maxIndex[0U]) {
            yR_1d = table[offset_0d];
          } else {
            yR_1d = (table[offset_0d + 1U] - table[offset_0d]) * frac[0U] +
              table[offset_0d];
          }

          yL_1d_1 += (yR_1d - yL_1d_1) * frac[1U];
        }

        yL_1d_0 += (yL_1d_1 - yL_1d_0) * frac[2U];
      }

      yL_1d += (yL_1d_0 - yL_1d) * frac[3U];
    }

    y += (yL_1d - y) * frac[4U];
  }

  return y;
}

real_T intrp8d_s32dl_pw(const int32_T bpIndex[], const real_T frac[], const
  real_T table[], const uint32_T stride[])
{
  real_T yL;
  real_T yR;
  uint32_T offset;

  /* Interpolation 8-D
     Interpolation method: 'Linear'
     Use last breakpoint for index at or above upper limit: 'off'
     Overflow mode: 'portable wrapping'
   */
  offset = (uint32_T)(bpIndex[7U] * (int32_T)stride[7U]);
  yL = intrp7d_s32dl_pw(bpIndex, frac, &table[offset], stride);
  offset += stride[7U];
  yR = intrp7d_s32dl_pw(bpIndex, frac, &table[offset], stride);
  return (yR - yL) * frac[7U] + yL;
}

real_T intrp7d_s32dl_pw(const int32_T bpIndex[], const real_T frac[], const
  real_T table[], const uint32_T stride[])
{
  real_T yL;
  real_T yR;
  uint32_T offset;

  /* Interpolation 7-D
     Interpolation method: 'Linear'
     Use last breakpoint for index at or above upper limit: 'off'
     Overflow mode: 'portable wrapping'
   */
  offset = (uint32_T)(bpIndex[6U] * (int32_T)stride[6U]);
  yL = intrp6d_s32dl_pw(bpIndex, frac, &table[offset], stride);
  offset += stride[6U];
  yR = intrp6d_s32dl_pw(bpIndex, frac, &table[offset], stride);
  return (yR - yL) * frac[6U] + yL;
}

real_T intrp6d_s32dl_pw(const int32_T bpIndex[], const real_T frac[], const
  real_T table[], const uint32_T stride[])
{
  real_T yL;
  real_T yR;
  uint32_T offset;

  /* Interpolation 6-D
     Interpolation method: 'Linear'
     Use last breakpoint for index at or above upper limit: 'off'
     Overflow mode: 'portable wrapping'
   */
  offset = (uint32_T)(bpIndex[5U] * (int32_T)stride[5U]);
  yL = intrp5d_s32dl_pw(bpIndex, frac, &table[offset], stride);
  offset += stride[5U];
  yR = intrp5d_s32dl_pw(bpIndex, frac, &table[offset], stride);
  return (yR - yL) * frac[5U] + yL;
}

real_T intrp5d_s32dl_pw(const int32_T bpIndex[], const real_T frac[], const
  real_T table[], const uint32_T stride[])
{
  real_T yL_4d;
  uint32_T offset_4d;
  real_T yL_3d;
  real_T yL_2d;
  real_T yL_1d;
  uint32_T offset_0d;
  uint32_T offset_1d;

  /* Interpolation 5-D
     Interpolation method: 'Linear'
     Use last breakpoint for index at or above upper limit: 'off'
     Overflow mode: 'portable wrapping'
   */
  offset_4d = ((((uint32_T)(bpIndex[4U] * (int32_T)stride[4U]) + bpIndex[3U] *
                 (int32_T)stride[3U]) + bpIndex[2U] * (int32_T)stride[2U]) +
               bpIndex[1U] * (int32_T)stride[1U]) + bpIndex[0U];
  yL_1d = (table[offset_4d + 1U] - table[offset_4d]) * frac[0U] +
    table[offset_4d];
  offset_0d = offset_4d + stride[1U];
  yL_2d = (((table[offset_0d + 1U] - table[offset_0d]) * frac[0U] +
            table[offset_0d]) - yL_1d) * frac[1U] + yL_1d;
  offset_1d = offset_4d + stride[2U];
  yL_1d = (table[offset_1d + 1U] - table[offset_1d]) * frac[0U] +
    table[offset_1d];
  offset_0d = offset_1d + stride[1U];
  yL_3d = (((((table[offset_0d + 1U] - table[offset_0d]) * frac[0U] +
              table[offset_0d]) - yL_1d) * frac[1U] + yL_1d) - yL_2d) * frac[2U]
    + yL_2d;
  offset_1d = offset_4d + stride[3U];
  yL_1d = (table[offset_1d + 1U] - table[offset_1d]) * frac[0U] +
    table[offset_1d];
  offset_0d = offset_1d + stride[1U];
  yL_2d = (((table[offset_0d + 1U] - table[offset_0d]) * frac[0U] +
            table[offset_0d]) - yL_1d) * frac[1U] + yL_1d;
  offset_1d += stride[2U];
  yL_1d = (table[offset_1d + 1U] - table[offset_1d]) * frac[0U] +
    table[offset_1d];
  offset_0d = offset_1d + stride[1U];
  yL_4d = (((((((table[offset_0d + 1U] - table[offset_0d]) * frac[0U] +
                table[offset_0d]) - yL_1d) * frac[1U] + yL_1d) - yL_2d) * frac
            [2U] + yL_2d) - yL_3d) * frac[3U] + yL_3d;
  offset_4d += stride[4U];
  yL_1d = (table[offset_4d + 1U] - table[offset_4d]) * frac[0U] +
    table[offset_4d];
  offset_0d = offset_4d + stride[1U];
  yL_2d = (((table[offset_0d + 1U] - table[offset_0d]) * frac[0U] +
            table[offset_0d]) - yL_1d) * frac[1U] + yL_1d;
  offset_1d = offset_4d + stride[2U];
  yL_1d = (table[offset_1d + 1U] - table[offset_1d]) * frac[0U] +
    table[offset_1d];
  offset_0d = offset_1d + stride[1U];
  yL_3d = (((((table[offset_0d + 1U] - table[offset_0d]) * frac[0U] +
              table[offset_0d]) - yL_1d) * frac[1U] + yL_1d) - yL_2d) * frac[2U]
    + yL_2d;
  offset_1d = offset_4d + stride[3U];
  yL_1d = (table[offset_1d + 1U] - table[offset_1d]) * frac[0U] +
    table[offset_1d];
  offset_0d = offset_1d + stride[1U];
  yL_2d = (((table[offset_0d + 1U] - table[offset_0d]) * frac[0U] +
            table[offset_0d]) - yL_1d) * frac[1U] + yL_1d;
  offset_1d += stride[2U];
  yL_1d = (table[offset_1d + 1U] - table[offset_1d]) * frac[0U] +
    table[offset_1d];
  offset_0d = offset_1d + stride[1U];
  return (((((((((table[offset_0d + 1U] - table[offset_0d]) * frac[0U] +
                 table[offset_0d]) - yL_1d) * frac[1U] + yL_1d) - yL_2d) * frac
             [2U] + yL_2d) - yL_3d) * frac[3U] + yL_3d) - yL_4d) * frac[4U] +
    yL_4d;
}

uint32_T plook_evenca(real_T u, real_T bp0, real_T bpSpace, uint32_T maxIndex,
                      real_T *fraction)
{
  uint32_T bpIndex;
  real_T invSpc;
  real_T fbpIndex;

  /* Prelookup - Index and Fraction
     Index Search method: 'even'
     Extrapolation method: 'Clip'
     Use previous index: 'off'
     Use last breakpoint for index at or above upper limit: 'on'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u <= bp0) {
    bpIndex = 0U;
    *fraction = 0.0;
  } else {
    invSpc = 1.0 / bpSpace;
    fbpIndex = (u - bp0) * invSpc;
    if (fbpIndex < maxIndex) {
      bpIndex = (uint32_T)fbpIndex;
      *fraction = (u - ((real_T)bpIndex * bpSpace + bp0)) * invSpc;
    } else {
      bpIndex = maxIndex;
      *fraction = 0.0;
    }
  }

  return bpIndex;
}

real_T intrp3d_la_pw(const uint32_T bpIndex[], const real_T frac[], const real_T
                     table[], const uint32_T stride[], const uint32_T maxIndex[])
{
  real_T y;
  uint32_T offset_2d;
  real_T yR_1d;
  uint32_T offset_0d;
  real_T yL_1d;

  /* Interpolation 3-D
     Interpolation method: 'Linear'
     Use last breakpoint for index at or above upper limit: 'on'
     Overflow mode: 'portable wrapping'
   */
  offset_2d = (bpIndex[2U] * stride[2U] + bpIndex[1U] * stride[1U]) + bpIndex[0U];
  if (bpIndex[0U] == maxIndex[0U]) {
    y = table[offset_2d];
  } else {
    y = (table[offset_2d + 1U] - table[offset_2d]) * frac[0U] + table[offset_2d];
  }

  if (bpIndex[1U] == maxIndex[1U]) {
  } else {
    offset_0d = offset_2d + stride[1U];
    if (bpIndex[0U] == maxIndex[0U]) {
      yR_1d = table[offset_0d];
    } else {
      yR_1d = (table[offset_0d + 1U] - table[offset_0d]) * frac[0U] +
        table[offset_0d];
    }

    y += (yR_1d - y) * frac[1U];
  }

  if (bpIndex[2U] == maxIndex[2U]) {
  } else {
    offset_2d += stride[2U];
    if (bpIndex[0U] == maxIndex[0U]) {
      yL_1d = table[offset_2d];
    } else {
      yL_1d = (table[offset_2d + 1U] - table[offset_2d]) * frac[0U] +
        table[offset_2d];
    }

    if (bpIndex[1U] == maxIndex[1U]) {
    } else {
      offset_0d = offset_2d + stride[1U];
      if (bpIndex[0U] == maxIndex[0U]) {
        yR_1d = table[offset_0d];
      } else {
        yR_1d = (table[offset_0d + 1U] - table[offset_0d]) * frac[0U] +
          table[offset_0d];
      }

      yL_1d += (yR_1d - yL_1d) * frac[1U];
    }

    y += (yL_1d - y) * frac[2U];
  }

  return y;
}

uint32_T plook_binca(real_T u, const real_T bp[], uint32_T maxIndex, real_T
                     *fraction)
{
  uint32_T bpIndex;

  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Clip'
     Use previous index: 'off'
     Use last breakpoint for index at or above upper limit: 'on'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u <= bp[0U]) {
    bpIndex = 0U;
    *fraction = 0.0;
  } else if (u < bp[maxIndex]) {
    bpIndex = binsearch_u32d(u, bp, maxIndex >> 1U, maxIndex);
    *fraction = (u - bp[bpIndex]) / (bp[bpIndex + 1U] - bp[bpIndex]);
  } else {
    bpIndex = maxIndex;
    *fraction = 0.0;
  }

  return bpIndex;
}

int32_T binsearch_s32d_prevIdx(real_T u, const real_T bp[], uint32_T startIndex,
  uint32_T maxIndex)
{
  uint32_T iRght;
  uint32_T iLeft;
  uint32_T bpIdx;
  uint32_T found;

  /* Binary Search using Previous Index */
  bpIdx = startIndex;
  iLeft = 0U;
  iRght = maxIndex;
  found = 0U;
  while (found == 0U) {
    if (u < bp[bpIdx]) {
      iRght = bpIdx - 1U;
      bpIdx = (iRght + iLeft) >> 1U;
    } else if (u < bp[bpIdx + 1U]) {
      found = 1U;
    } else {
      iLeft = bpIdx + 1U;
      bpIdx = (iRght + iLeft) >> 1U;
    }
  }

  return (int32_T)bpIdx;
}

uint32_T binsearch_u32d_prevIdx(real_T u, const real_T bp[], uint32_T startIndex,
  uint32_T maxIndex)
{
  uint32_T bpIndex;
  uint32_T iRght;
  uint32_T iLeft;
  uint32_T found;

  /* Binary Search using Previous Index */
  bpIndex = startIndex;
  iLeft = 0U;
  iRght = maxIndex;
  found = 0U;
  while (found == 0U) {
    if (u < bp[bpIndex]) {
      iRght = bpIndex - 1U;
      bpIndex = (iRght + iLeft) >> 1U;
    } else if (u < bp[bpIndex + 1U]) {
      found = 1U;
    } else {
      iLeft = bpIndex + 1U;
      bpIndex = (iRght + iLeft) >> 1U;
    }
  }

  return bpIndex;
}

uint32_T binsearch_u32d(real_T u, const real_T bp[], uint32_T startIndex,
  uint32_T maxIndex)
{
  uint32_T bpIndex;
  uint32_T iRght;
  uint32_T bpIdx;

  /* Binary Search */
  bpIdx = startIndex;
  bpIndex = 0U;
  iRght = maxIndex;
  while (iRght - bpIndex > 1U) {
    if (u < bp[bpIdx]) {
      iRght = bpIdx;
    } else {
      bpIndex = bpIdx;
    }

    bpIdx = (iRght + bpIndex) >> 1U;
  }

  return bpIndex;
}

/*
 *   This function updates active task flag for each subrate.
 * The function is called at model base rate, hence the
 * generated code self-manages all its subrates.
 */
static void rate_scheduler(void)
{
  /* Compute which subrates run during the next base time step.  Subrates
   * are an integer multiple of the base rate counter.  Therefore, the subtask
   * counter is reset when it reaches its limit (zero means run).
   */
  (mainV03_56_M->Timing.TaskCounters.TID[2])++;
  if ((mainV03_56_M->Timing.TaskCounters.TID[2]) > 4) {/* Sample time: [0.1s, 0.0s] */
    mainV03_56_M->Timing.TaskCounters.TID[2] = 0;
  }

  mainV03_56_M->Timing.sampleHits[2] = (mainV03_56_M->Timing.TaskCounters.TID[2]
    == 0);
}

/* Simplified version of numjac.cpp, for use with RTW. */
void local_numjac( RTWSolverInfo *si, real_T *y, const real_T *Fty, real_T *fac,
                  real_T *dFdy )
{
  /* constants */
  real_T THRESH = 1e-6;
  real_T EPS = 2.2e-16;                /* utGetEps(); */
  real_T BL = pow(EPS, 0.75);
  real_T BU = pow(EPS, 0.25);
  real_T FACMIN = pow(EPS, 0.78);
  real_T FACMAX = 0.1;
  int_T nx = 140;
  real_T *x = rtsiGetContStates(si);
  real_T del;
  real_T difmax;
  real_T FdelRowmax;
  real_T temp;
  real_T Fdiff;
  real_T maybe;
  real_T xscale;
  real_T fscale;
  real_T *p;
  int_T rowmax;
  int_T i,j;
  if (x != y)
    (void) memcpy(x, y,
                  (uint_T)nx*sizeof(real_T));
  rtsiSetSolverComputingJacobian(si,true);
  for (p = dFdy, j = 0; j < nx; j++, p += nx) {
    /* Select an increment del for a difference approximation to
       column j of dFdy.  The vector fac accounts for experience
       gained in previous calls to numjac. */
    xscale = fabs(x[j]);
    if (xscale < THRESH)
      xscale = THRESH;
    temp = (x[j] + fac[j]*xscale);
    del = temp - y[j];
    while (del == 0.0) {
      if (fac[j] < FACMAX) {
        fac[j] *= 100.0;
        if (fac[j] > FACMAX)
          fac[j] = FACMAX;
        temp = (x[j] + fac[j]*xscale);
        del = temp - x[j];
      } else {
        del = THRESH;                  /* thresh is nonzero */
        break;
      }
    }

    /* Keep del pointing into region. */
    if (Fty[j] >= 0.0)
      del = fabs(del);
    else
      del = -fabs(del);

    /* Form a difference approximation to column j of dFdy. */
    temp = x[j];
    x[j] += del;
    mainV03_56_output();
    rtsiSetdX(si,p);
    mainV03_56_derivatives();
    x[j] = temp;
    difmax = 0.0;
    rowmax = 0;
    FdelRowmax = p[0];
    temp = 1.0 / del;
    for (i = 0; i < nx; i++) {
      Fdiff = p[i] - Fty[i];
      maybe = fabs(Fdiff);
      if (maybe > difmax) {
        difmax = maybe;
        rowmax = i;
        FdelRowmax = p[i];
      }

      p[i] = temp * Fdiff;
    }

    /* Adjust fac for next call to numjac. */
    if (((FdelRowmax != 0.0) && (Fty[rowmax] != 0.0)) || (difmax == 0.0)) {
      fscale = fabs(FdelRowmax);
      if (fscale < fabs(Fty[rowmax]))
        fscale = fabs(Fty[rowmax]);
      if (difmax <= BL*fscale) {
        /* The difference is small, so increase the increment. */
        fac[j] *= 10.0;
        if (fac[j] > FACMAX)
          fac[j] = FACMAX;
      } else if (difmax > BU*fscale) {
        /* The difference is large, so reduce the increment. */
        fac[j] *= 0.1;
        if (fac[j] < FACMIN)
          fac[j] = FACMIN;
      }
    }
  }

  rtsiSetSolverComputingJacobian(si,false);
}                                      /* end local_numjac */

/*
 * This function updates continuous states using the ODE14x fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  /* Solver Matrices */
  static int_T rt_ODE14x_N[4] = { 12, 8, 6, 4 };

  time_T t0 = rtsiGetT(si);
  time_T t1 = t0;
  time_T h = rtsiGetStepSize(si);
  real_T *x1 = rtsiGetContStates(si);
  int_T order = rtsiGetSolverExtrapolationOrder(si);
  int_T numIter = rtsiGetSolverNumberNewtonIterations(si);
  ODE14X_IntgData *id = (ODE14X_IntgData *)rtsiGetSolverData(si);
  real_T *x0 = id->x0;
  real_T *f0 = id->f0;
  real_T *x1start = id->x1start;
  real_T *f1 = id->f1;
  real_T *Delta = id->Delta;
  real_T *E = id->E;
  real_T *fac = id->fac;
  real_T *dfdx = id->DFDX;
  real_T *W = id->W;
  int_T *pivots = id->pivots;
  real_T *xtmp = id->xtmp;
  real_T *ztmp = id->ztmp;
  int_T *N = &(rt_ODE14x_N[0]);
  int_T i,j,k,iter;
  int_T nx = 140;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) memcpy(x0, x1,
                (uint_T)nx*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */

  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  mainV03_56_derivatives();
  local_numjac(si,x0,f0,fac,dfdx );
  for (j = 0; j < order; j++) {
    real_T *p;
    real_T hN = h/N[j];

    /* Get the iteration matrix and solution at t0 */

    /* [L,U] = lu(M - hN*J) */
    (void) memcpy(W, dfdx,
                  (uint_T)nx*nx*sizeof(real_T));
    for (p = W, i = 0; i < nx*nx; i++, p++) {
      *p *= (-hN);
    }

    for (p = W, i = 0; i < nx; i++, p += (nx+1)) {
      *p += 1.0;
    }

    rt_lu_real(W, nx,
               pivots);

    /* First Newton's iteration at t0. */
    /* rhs = hN*f0 */
    for (i = 0; i < nx; i++) {
      Delta[i] = hN*f0[i];
    }

    /* Delta = (U \ (L \ rhs)) */
    rt_ForwardSubstitutionRR_Dbl(W, Delta,
      f1, nx,
      1, pivots,
      1);
    rt_BackwardSubstitutionRR_Dbl(W+nx*nx-1, f1+nx-1,
      Delta, nx,
      1, 0);

    /* ytmp = y0 + Delta
       ztmp = (ytmp-y0)/h
     */
    (void) memcpy(x1, x0,
                  (uint_T)nx*sizeof(real_T));
    for (i = 0; i < nx; i++) {
      x1[i] += Delta[i];
      ztmp[i] = Delta[i]/hN;
    }

    /* Additional Newton's iterations, if desired.
       for iter = 2:NewtIter
       rhs = hN*feval(odefun,tn,ytmp,extraArgs{:}) - M*(ytmp - yn);
       if statedepM   % only for state dep. Mdel ~= 0
       Mdel = M - feval(massfun,tn,ytmp);
       rhs = rhs + Mdel*ztmp*h;
       end
       Delta = ( U \ ( L \ rhs ) );
       ytmp = ytmp + Delta;
       ztmp = (ytmp - yn)/h
       end
     */
    rtsiSetT(si, t0);
    rtsiSetdX(si, f1);
    for (iter = 1; iter < numIter; iter++) {
      mainV03_56_output();
      mainV03_56_derivatives();
      for (i = 0; i < nx; i++) {
        Delta[i] = hN*f1[i];
        xtmp[i] = x1[i] - x0[i];
      }

      /* rhs = hN*f(tn,ytmp) - (ytmp-yn) */
      for (i = 0; i < nx; i++) {
        Delta[i] -= xtmp[i];
      }

      rt_ForwardSubstitutionRR_Dbl(W, Delta,
        f1, nx,
        1, pivots,
        1);
      rt_BackwardSubstitutionRR_Dbl(W+nx*nx-1, f1+nx-1,
        Delta, nx,
        1, 0);

      /* ytmp = ytmp + delta
         ztmp = (ytmp - yn)/h
       */
      for (i = 0; i < nx; i++) {
        x1[i] += Delta[i];
        ztmp[i] = (x1[i] - x0[i])/hN;
      }
    }

    /* Steps from t0+hN to t1 -- subintegration of N(j) steps for extrapolation
       ttmp = t0;
       for i = 2:N(j)
       ttmp = ttmp + hN
       ytmp0 = ytmp;
       for iter = 1:NewtIter
       rhs = (ytmp0 - ytmp) + hN*feval(odefun,ttmp,ytmp,extraArgs{:});
       Delta = ( U \ ( L \ rhs ) );
       ytmp = ytmp + Delta;
       end
       end
     */
    for (k = 1; k < N[j]; k++) {
      t1 = t0 + k*hN;
      (void) memcpy(x1start, x1,
                    (uint_T)nx*sizeof(real_T));
      rtsiSetT(si, t1);
      rtsiSetdX(si, f1);
      for (iter = 0; iter < numIter; iter++) {
        mainV03_56_output();
        mainV03_56_derivatives();
        if (iter == 0) {
          for (i = 0; i < nx; i++) {
            Delta[i] = hN*f1[i];
          }
        } else {
          for (i = 0; i < nx; i++) {
            Delta[i] = hN*f1[i];
            xtmp[i] = (x1[i]-x1start[i]);
          }

          /* rhs = hN*f(tn,ytmp) - M*(ytmp-yn) */
          for (i = 0; i < nx; i++) {
            Delta[i] -= xtmp[i];
          }
        }

        rt_ForwardSubstitutionRR_Dbl(W, Delta,
          f1, nx,
          1, pivots,
          1);
        rt_BackwardSubstitutionRR_Dbl(W+nx*nx-1, f1+nx-1,
          Delta, nx,
          1, 0);

        /* ytmp = ytmp + Delta
           ztmp = (ytmp - ytmp0)/h
         */
        for (i = 0; i < nx; i++) {
          x1[i] += Delta[i];
          ztmp[i] = (x1[i] - x1start[i])/hN;
        }
      }
    }

    /* Extrapolate to order j
       E(:,j) = ytmp
       for k = j:-1:2
       coef = N(k-1)/(N(j) - N(k-1))
       E(:,k-1) = E(:,k) + coef*( E(:,k) - E(:,k-1) )
       end
     */
    (void) memcpy(&(E[nx*j]), x1,
                  (uint_T)nx*sizeof(real_T));
    for (k = j; k > 0; k--) {
      real_T coef = (real_T)(N[k-1]) / (N[j]-N[k-1]);
      for (i = 0; i < nx; i++) {
        x1[i] = E[nx*k+i] + coef*(E[nx*k+i] - E[nx*(k-1)+i]);
      }

      (void) memcpy(&(E[nx*(k-1)]), x1,
                    (uint_T)nx*sizeof(real_T));
    }
  }

  /* x1 = E(:,1); */
  (void) memcpy(x1, E,
                (uint_T)nx*sizeof(real_T));

  /* t1 = t0 + h; */
  rtsiSetT(si,rtsiGetSolverStopTime(si));
  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/*
 * System initialize for enable system:
 *    '<S47>/Distance into gust (y)'
 *    '<S47>/Distance into gust (z)'
 */
void mainV03_56_Distanceintogusty_Init(B_Distanceintogusty_mainV03_56_T *localB,
  P_Distanceintogusty_mainV03_56_T *localP, X_Distanceintogusty_mainV03_56_T
  *localX)
{
  /* InitializeConditions for Integrator: '<S51>/Distance into Gust (x) (Limited to gust length d) ' */
  localX->DistanceintoGustxLimitedtogustlengthd_CSTATE =
    localP->DistanceintoGustxLimitedtogustlengthd_IC;

  /* SystemInitialize for Outport: '<S51>/x' */
  localB->DistanceintoGustxLimitedtogustlengthd = localP->x_Y0;
}

/*
 * System reset for enable system:
 *    '<S47>/Distance into gust (y)'
 *    '<S47>/Distance into gust (z)'
 */
void mainV03_56_Distanceintogusty_Reset(P_Distanceintogusty_mainV03_56_T *localP,
  X_Distanceintogusty_mainV03_56_T *localX)
{
  /* InitializeConditions for Integrator: '<S51>/Distance into Gust (x) (Limited to gust length d) ' */
  localX->DistanceintoGustxLimitedtogustlengthd_CSTATE =
    localP->DistanceintoGustxLimitedtogustlengthd_IC;
}

/*
 * Disable for enable system:
 *    '<S47>/Distance into gust (y)'
 *    '<S47>/Distance into gust (z)'
 */
void mainV03_56_Distanceintogusty_Disable(DW_Distanceintogusty_mainV03_56_T
  *localDW)
{
  localDW->Distanceintogusty_MODE = false;
}

/*
 * Outputs for enable system:
 *    '<S47>/Distance into gust (y)'
 *    '<S47>/Distance into gust (z)'
 */
void mainV03_56_Distanceintogusty(RT_MODEL_mainV03_56_T * const mainV03_56_M,
  boolean_T rtu_Enable, B_Distanceintogusty_mainV03_56_T *localB,
  DW_Distanceintogusty_mainV03_56_T *localDW, P_Distanceintogusty_mainV03_56_T
  *localP, X_Distanceintogusty_mainV03_56_T *localX, real_T rtp_d_m)
{
  /* Outputs for Enabled SubSystem: '<S47>/Distance into gust (y)' incorporates:
   *  EnablePort: '<S51>/Enable'
   */
  if ((rtmIsMajorTimeStep(mainV03_56_M) &&
       mainV03_56_M->Timing.TaskCounters.TID[1] == 0) && rtmIsMajorTimeStep
      (mainV03_56_M)) {
    if (rtu_Enable) {
      if (!localDW->Distanceintogusty_MODE) {
        mainV03_56_Distanceintogusty_Reset(localP, localX);
        localDW->Distanceintogusty_MODE = true;
      }
    } else {
      if (localDW->Distanceintogusty_MODE) {
        mainV03_56_Distanceintogusty_Disable(localDW);
      }
    }
  }

  if (localDW->Distanceintogusty_MODE) {
    /* Integrator: '<S51>/Distance into Gust (x) (Limited to gust length d) ' */
    /* Limited  Integrator  */
    if (localX->DistanceintoGustxLimitedtogustlengthd_CSTATE >= rtp_d_m) {
      localX->DistanceintoGustxLimitedtogustlengthd_CSTATE = rtp_d_m;
    } else {
      if (localX->DistanceintoGustxLimitedtogustlengthd_CSTATE <=
          localP->DistanceintoGustxLimitedtogustlengthd_LowerSat) {
        localX->DistanceintoGustxLimitedtogustlengthd_CSTATE =
          localP->DistanceintoGustxLimitedtogustlengthd_LowerSat;
      }
    }

    localB->DistanceintoGustxLimitedtogustlengthd =
      localX->DistanceintoGustxLimitedtogustlengthd_CSTATE;

    /* End of Integrator: '<S51>/Distance into Gust (x) (Limited to gust length d) ' */
  }

  /* End of Outputs for SubSystem: '<S47>/Distance into gust (y)' */
}

/*
 * Derivatives for enable system:
 *    '<S47>/Distance into gust (y)'
 *    '<S47>/Distance into gust (z)'
 */
void mainV03_56_Distanceintogusty_Deriv(real_T rtu_V,
  DW_Distanceintogusty_mainV03_56_T *localDW, P_Distanceintogusty_mainV03_56_T
  *localP, X_Distanceintogusty_mainV03_56_T *localX,
  XDot_Distanceintogusty_mainV03_56_T *localXdot, real_T rtp_d_m)
{
  boolean_T lsat;
  boolean_T usat;
  if (localDW->Distanceintogusty_MODE) {
    /* Derivatives for Integrator: '<S51>/Distance into Gust (x) (Limited to gust length d) ' */
    lsat = (localX->DistanceintoGustxLimitedtogustlengthd_CSTATE <=
            localP->DistanceintoGustxLimitedtogustlengthd_LowerSat);
    usat = (localX->DistanceintoGustxLimitedtogustlengthd_CSTATE >= rtp_d_m);
    if (((!lsat) && (!usat)) || (lsat && (rtu_V > 0.0)) || (usat && (rtu_V < 0.0)))
    {
      localXdot->DistanceintoGustxLimitedtogustlengthd_CSTATE = rtu_V;
    } else {
      /* in saturation */
      localXdot->DistanceintoGustxLimitedtogustlengthd_CSTATE = 0.0;
    }

    /* End of Derivatives for Integrator: '<S51>/Distance into Gust (x) (Limited to gust length d) ' */
  } else {
    localXdot->DistanceintoGustxLimitedtogustlengthd_CSTATE = 0.0;
  }
}

real_T rt_urand_Upu32_Yd_f_pw_snf(uint32_T *u)
{
  uint32_T lo;
  uint32_T hi;

  /* Uniform random number generator (random number between 0 and 1)

     #define IA      16807                      magic multiplier = 7^5
     #define IM      2147483647                 modulus = 2^31-1
     #define IQ      127773                     IM div IA
     #define IR      2836                       IM modulo IA
     #define S       4.656612875245797e-10      reciprocal of 2^31-1
     test = IA * (seed % IQ) - IR * (seed/IQ)
     seed = test < 0 ? (test + IM) : test
     return (seed*S)
   */
  lo = *u % 127773U * 16807U;
  hi = *u / 127773U * 2836U;
  if (lo < hi) {
    *u = 2147483647U - (hi - lo);
  } else {
    *u = lo - hi;
  }

  return (real_T)*u * 4.6566128752457969E-10;
}

real_T rt_nrand_Upu32_Yd_f_pw_snf(uint32_T *u)
{
  real_T y;
  real_T sr;
  real_T si;

  /* Normal (Gaussian) random number generator */
  do {
    sr = 2.0 * rt_urand_Upu32_Yd_f_pw_snf(u) - 1.0;
    si = 2.0 * rt_urand_Upu32_Yd_f_pw_snf(u) - 1.0;
    si = sr * sr + si * si;
  } while (si > 1.0);

  y = sqrt(-2.0 * log(si) / si) * sr;
  return y;
}

real_T rt_roundd_snf(real_T u)
{
  real_T y;
  if (fabs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = floor(u + 0.5);
    } else if (u > -0.5) {
      y = u * 0.0;
    } else {
      y = ceil(u - 0.5);
    }
  } else {
    y = u;
  }

  return y;
}

real_T rt_modd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T tmp;
  if (u1 == 0.0) {
    y = u0;
  } else if (!((!rtIsNaN(u0)) && (!rtIsInf(u0)) && ((!rtIsNaN(u1)) && (!rtIsInf
                (u1))))) {
    y = (rtNaN);
  } else {
    tmp = u0 / u1;
    if (u1 <= floor(u1)) {
      y = u0 - floor(tmp) * u1;
    } else if (fabs(tmp - rt_roundd_snf(tmp)) <= DBL_EPSILON * fabs(tmp)) {
      y = 0.0;
    } else {
      y = (tmp - floor(tmp)) * u1;
    }
  }

  return y;
}

real_T rt_atan2d_snf(real_T u0, real_T u1)
{
  real_T y;
  int32_T u0_0;
  int32_T u1_0;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    if (u0 > 0.0) {
      u0_0 = 1;
    } else {
      u0_0 = -1;
    }

    if (u1 > 0.0) {
      u1_0 = 1;
    } else {
      u1_0 = -1;
    }

    y = atan2(u0_0, u1_0);
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = atan2(u0, u1);
  }

  return y;
}

real_T rt_powd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T tmp;
  real_T tmp_0;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else {
    tmp = fabs(u0);
    tmp_0 = fabs(u1);
    if (rtIsInf(u1)) {
      if (tmp == 1.0) {
        y = (rtNaN);
      } else if (tmp > 1.0) {
        if (u1 > 0.0) {
          y = (rtInf);
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = (rtInf);
      }
    } else if (tmp_0 == 0.0) {
      y = 1.0;
    } else if (tmp_0 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = (rtNaN);
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

void rt_mrdivide_U1d1x3_U2d3x3_Yd1x3_snf(const real_T u0[3], const real_T u1[9],
  real_T y[3])
{
  real_T A[9];
  int32_T r1;
  int32_T r2;
  int32_T r3;
  real_T maxval;
  real_T a21;
  int32_T rtemp;
  memcpy(&A[0], &u1[0], 9U * sizeof(real_T));
  r1 = 0;
  r2 = 1;
  r3 = 2;
  maxval = fabs(u1[0]);
  a21 = fabs(u1[1]);
  if (a21 > maxval) {
    maxval = a21;
    r1 = 1;
    r2 = 0;
  }

  if (fabs(u1[2]) > maxval) {
    r1 = 2;
    r2 = 1;
    r3 = 0;
  }

  A[r2] = u1[r2] / u1[r1];
  A[r3] /= A[r1];
  A[3 + r2] -= A[3 + r1] * A[r2];
  A[3 + r3] -= A[3 + r1] * A[r3];
  A[6 + r2] -= A[6 + r1] * A[r2];
  A[6 + r3] -= A[6 + r1] * A[r3];
  if (fabs(A[3 + r3]) > fabs(A[3 + r2])) {
    rtemp = r2 + 1;
    r2 = r3;
    r3 = rtemp - 1;
  }

  A[3 + r3] /= A[3 + r2];
  A[6 + r3] -= A[3 + r3] * A[6 + r2];
  y[r1] = u0[0] / A[r1];
  y[r2] = u0[1] - A[3 + r1] * y[r1];
  y[r3] = u0[2] - A[6 + r1] * y[r1];
  y[r2] /= A[3 + r2];
  y[r3] -= A[6 + r2] * y[r2];
  y[r3] /= A[6 + r3];
  y[r2] -= A[3 + r3] * y[r3];
  y[r1] -= y[r3] * A[r3];
  y[r1] -= y[r2] * A[r2];
}

/* Model output function */
void mainV03_56_output(void)
{
  /* local block i/o variables */
  real_T rtb_FromWs[6];
  real_T rtb_sincos_o2_p[3];
  real_T rtb_UnaryMinus_o;
  real_T rtb_UnaryMinus_j;
  boolean_T rtb_Compare_a;

  /* local scratch DWork variables */
  int32_T ForEach_itr;
  real_T *lastU;
  real_T riseValLimit;
  real_T (*lastU_0)[3];
  real_T rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[2];
  int32_T rtb_PreLookUpIndexSearch_o1;
  real_T rtb_Rn;
  real_T rtb_Sum1_no;
  real_T rtb_Switch;
  real_T rtb_Abs1;
  boolean_T rtb_Compare_e;
  int8_T rtAction;
  int32_T s106_iter;
  real_T rtb_WhiteNoise_o[3];
  real_T rtb_VectorConcatenate_i[9];
  real_T rtb_VectorConcatenate_n[9];
  int32_T s155_iter;
  real_T rtb_Assignment1[11];
  real_T rtb_Assignment[11];
  int32_T s114_iter;
  real_T rtb_BusCreator3_n_Throttle2;
  real_T rtb_BusCreator3_n_Throttle3;
  real_T rtb_BusCreator3_n_Throttle4;
  real_T rtb_BusCreator3_n_Throttle5;
  real_T rtb_BusCreator1_Throttle1;
  real_T rtb_BusCreator1_Throttle2;
  real_T rtb_BusCreator1_Throttle3;
  real_T rtb_BusCreator1_Throttle4;
  real_T rtb_BusCreator1_Throttle5;
  real_T rtb_BusCreator3_o_Throttle1;
  real_T rtb_BusCreator3_o_Throttle2;
  real_T rtb_BusCreator3_o_Throttle3;
  real_T rtb_BusCreator3_o_Throttle4;
  real_T rtb_BusCreator3_o_Throttle5;
  int32_T rtb_Sum1_h;
  real_T rtb_Shiftright[3];
  real_T frac[2];
  uint32_T bpIndex[2];
  real_T frac_0[8];
  int32_T bpIndex_0[8];
  real_T frac_1[8];
  int32_T bpIndex_1[8];
  real_T frac_2[8];
  int32_T bpIndex_2[8];
  real_T frac_3[8];
  int32_T bpIndex_3[8];
  real_T frac_4[8];
  int32_T bpIndex_4[8];
  real_T frac_5[8];
  int32_T bpIndex_5[8];
  real_T frac_6[8];
  int32_T bpIndex_6[8];
  real_T frac_7[8];
  int32_T bpIndex_7[8];
  real_T frac_8[8];
  int32_T bpIndex_8[8];
  real_T frac_9[8];
  int32_T bpIndex_9[8];
  real_T frac_a[8];
  int32_T bpIndex_a[8];
  real_T frac_b[8];
  int32_T bpIndex_b[8];
  real_T frac_c[8];
  int32_T bpIndex_c[8];
  real_T frac_d[8];
  int32_T bpIndex_d[8];
  real_T frac_e[8];
  int32_T bpIndex_e[8];
  real_T frac_f[8];
  int32_T bpIndex_f[8];
  real_T frac_g[8];
  int32_T bpIndex_g[8];
  real_T frac_h[8];
  int32_T bpIndex_h[8];
  real_T frac_i[8];
  int32_T bpIndex_i[8];
  real_T frac_j[8];
  int32_T bpIndex_j[8];
  real_T frac_k[8];
  int32_T bpIndex_k[8];
  real_T frac_l[8];
  int32_T bpIndex_l[8];
  real_T frac_m[8];
  int32_T bpIndex_m[8];
  real_T frac_n[8];
  int32_T bpIndex_n[8];
  real_T frac_o[8];
  int32_T bpIndex_o[8];
  real_T frac_p[8];
  int32_T bpIndex_p[8];
  real_T frac_q[8];
  int32_T bpIndex_q[8];
  real_T frac_r[8];
  int32_T bpIndex_r[8];
  real_T frac_s[8];
  int32_T bpIndex_s[8];
  real_T frac_t[8];
  int32_T bpIndex_t[8];
  real_T frac_u[8];
  int32_T bpIndex_u[8];
  real_T frac_v[8];
  int32_T bpIndex_v[8];
  real_T frac_w[8];
  int32_T bpIndex_w[8];
  real_T frac_x[8];
  int32_T bpIndex_x[8];
  real_T frac_y[8];
  int32_T bpIndex_y[8];
  real_T frac_z[8];
  int32_T bpIndex_z[8];
  real_T frac_10[8];
  int32_T bpIndex_10[8];
  real_T frac_11[8];
  int32_T bpIndex_11[8];
  real_T frac_12[8];
  int32_T bpIndex_12[8];
  real_T frac_13[3];
  uint32_T bpIndex_13[3];
  real_T frac_14[3];
  uint32_T bpIndex_14[3];
  real_T frac_15[3];
  uint32_T bpIndex_15[3];
  real_T frac_16[3];
  uint32_T bpIndex_16[3];
  real_T frac_17[3];
  uint32_T bpIndex_17[3];
  real_T frac_18[3];
  uint32_T bpIndex_18[3];
  real_T frac_19[3];
  uint32_T bpIndex_19[3];
  real_T frac_1a[3];
  uint32_T bpIndex_1a[3];
  real_T frac_1b[3];
  uint32_T bpIndex_1b[3];
  real_T frac_1c[3];
  uint32_T bpIndex_1c[3];
  real_T frac_1d[3];
  uint32_T bpIndex_1d[3];
  real_T frac_1e[3];
  uint32_T bpIndex_1e[3];
  real_T frac_1f[3];
  uint32_T bpIndex_1f[3];
  real_T frac_1g[3];
  uint32_T bpIndex_1g[3];
  real_T frac_1h[3];
  uint32_T bpIndex_1h[3];
  real_T frac_1i[2];
  uint32_T bpIndex_1i[2];
  real_T frac_1j[2];
  uint32_T bpIndex_1j[2];
  real_T frac_1k[2];
  uint32_T bpIndex_1k[2];
  real_T frac_1l[2];
  uint32_T bpIndex_1l[2];
  real_T frac_1m[2];
  uint32_T bpIndex_1m[2];
  real_T rtb_MatrixConcatenation[18];
  real_T rtb_tc_old[169];
  real_T Assignment[169];
  real_T Assignment2[169];
  real_T rtb_sincos_o1_n_idx_0;
  real_T rtb_sincos_o1_n_idx_1;
  real_T rtb_IC4_idx_0;
  real_T rtb_IC4_idx_1;
  real_T rtb_IC4_idx_2;
  real_T rtb_IC3_idx_1;
  real_T rtb_IC3_idx_0;
  real_T rtb_IC3_idx_2;
  real_T rtb_sincos_o2_i_idx_0;
  real_T rtb_sincos_o2_i_idx_1;
  real_T rtb_sincos_o2_i_idx_2;
  real_T rtb_sincos_o2_g_idx_1;
  real_T rtb_sincos_o2_g_idx_0;
  real_T rtb_sincos_o1_idx_2;
  real_T rtb_sincos_o1_idx_1;
  real_T rtb_sincos_o2_g_idx_2;
  real_T rtb_sincos_o1_idx_0;
  real_T rtb_ImpAsg_InsertedFor_F_at_inport_0_idx_0;
  real_T rtb_ImpAsg_InsertedFor_F_at_inport_0_idx_1;
  real_T rtb_ImpAsg_InsertedFor_F_at_inport_0_idx_2;
  int64_T tmp;
  int32_T qY;
  int32_T s106_iter_0;
  int32_T rtb_Sum1_eh;
  if (rtmIsMajorTimeStep(mainV03_56_M)) {
    /* set solver stop time */
    if (!(mainV03_56_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&mainV03_56_M->solverInfo,
                            ((mainV03_56_M->Timing.clockTickH0 + 1) *
        mainV03_56_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&mainV03_56_M->solverInfo,
                            ((mainV03_56_M->Timing.clockTick0 + 1) *
        mainV03_56_M->Timing.stepSize0 + mainV03_56_M->Timing.clockTickH0 *
        mainV03_56_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(mainV03_56_M)) {
    mainV03_56_M->Timing.t[0] = rtsiGetT(&mainV03_56_M->solverInfo);
  }

  /* Integrator: '<S243>/xe,ye,ze' */
  mainV03_56_B.xeyeze[0] = mainV03_56_X.xeyeze_CSTATE[0];
  mainV03_56_B.xeyeze[1] = mainV03_56_X.xeyeze_CSTATE[1];
  mainV03_56_B.xeyeze[2] = mainV03_56_X.xeyeze_CSTATE[2];

  /* InitialCondition: '<S4>/IC4' */
  if ((mainV03_56_DW.IC4_FirstOutputTime == (rtMinusInf)) ||
      (mainV03_56_DW.IC4_FirstOutputTime == mainV03_56_M->Timing.t[0])) {
    mainV03_56_DW.IC4_FirstOutputTime = mainV03_56_M->Timing.t[0];
    rtb_IC4_idx_0 = mainV03_56_P.IC4_Value_e[0];
    rtb_IC4_idx_1 = mainV03_56_P.IC4_Value_e[1];
    rtb_IC4_idx_2 = mainV03_56_P.IC4_Value_e[2];
  } else {
    rtb_IC4_idx_0 = mainV03_56_B.xeyeze[0];
    rtb_IC4_idx_1 = mainV03_56_B.xeyeze[1];
    rtb_IC4_idx_2 = mainV03_56_B.xeyeze[2];
  }

  /* End of InitialCondition: '<S4>/IC4' */
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* UnitConversion: '<S327>/Unit Conversion' incorporates:
     *  Constant: '<S313>/ref_pos'
     */
    /* Unit Conversion - from: deg to: rad
       Expression: output = (0.0174533*input) + (0) */
    rtb_Rn = 0.017453292519943295 * mainV03_56_P.FlatEarthtoLLA_psi;

    /* Trigonometry: '<S313>/SinCos' */
    mainV03_56_B.SinCos_o1 = sin(rtb_Rn);
    mainV03_56_B.SinCos_o2 = cos(rtb_Rn);

    /* Sum: '<S331>/Sum' incorporates:
     *  Constant: '<S331>/Constant'
     *  Constant: '<S331>/f'
     */
    rtb_Rn = mainV03_56_P.f_Value - mainV03_56_P.Constant_Value_jd;

    /* Sqrt: '<S332>/sqrt' incorporates:
     *  Constant: '<S332>/Constant'
     *  Product: '<S332>/Product1'
     *  Sum: '<S332>/Sum1'
     */
    rtb_Rn = sqrt(mainV03_56_P.Constant_Value_hu - rtb_Rn * rtb_Rn);

    /* Switch: '<S324>/Switch' incorporates:
     *  Abs: '<S324>/Abs'
     *  Bias: '<S324>/Bias'
     *  Bias: '<S324>/Bias1'
     *  Constant: '<S244>/initial_pos'
     *  Constant: '<S324>/Constant2'
     *  Constant: '<S325>/Constant'
     *  Math: '<S324>/Math Function1'
     *  RelationalOperator: '<S325>/Compare'
     */
    if (fabs(mainV03_56_P.FlatEarthtoLLA_LL0[0]) >
        mainV03_56_P.CompareToConstant_const) {
      rtb_Switch = rt_modd_snf(mainV03_56_P.FlatEarthtoLLA_LL0[0] +
        mainV03_56_P.Bias_Bias_g4, mainV03_56_P.Constant2_Value_b) +
        mainV03_56_P.Bias1_Bias_b;
    } else {
      rtb_Switch = mainV03_56_P.FlatEarthtoLLA_LL0[0];
    }

    /* End of Switch: '<S324>/Switch' */

    /* Abs: '<S321>/Abs1' */
    rtb_Abs1 = fabs(rtb_Switch);

    /* RelationalOperator: '<S323>/Compare' incorporates:
     *  Constant: '<S323>/Constant'
     */
    rtb_Compare_e = (rtb_Abs1 > mainV03_56_P.CompareToConstant_const_c);

    /* Switch: '<S321>/Switch' incorporates:
     *  Bias: '<S321>/Bias'
     *  Bias: '<S321>/Bias1'
     *  Gain: '<S321>/Gain'
     *  Product: '<S321>/Divide1'
     */
    if (rtb_Compare_e) {
      /* Signum: '<S321>/Sign1' */
      if (rtb_Switch < 0.0) {
        rtb_Switch = -1.0;
      } else if (rtb_Switch > 0.0) {
        rtb_Switch = 1.0;
      } else {
        if (rtb_Switch == 0.0) {
          rtb_Switch = 0.0;
        }
      }

      /* End of Signum: '<S321>/Sign1' */
      mainV03_56_B.Switch = ((rtb_Abs1 + mainV03_56_P.Bias_Bias_p) *
        mainV03_56_P.Gain_Gain_dr + mainV03_56_P.Bias1_Bias_o) * rtb_Switch;
    } else {
      mainV03_56_B.Switch = rtb_Switch;
    }

    /* End of Switch: '<S321>/Switch' */

    /* UnitConversion: '<S329>/Unit Conversion' */
    /* Unit Conversion - from: deg to: rad
       Expression: output = (0.0174533*input) + (0) */
    rtb_Sum1_no = 0.017453292519943295 * mainV03_56_B.Switch;

    /* Trigonometry: '<S330>/Trigonometric Function1' */
    rtb_Abs1 = sin(rtb_Sum1_no);

    /* Sum: '<S330>/Sum1' incorporates:
     *  Constant: '<S330>/Constant'
     *  Product: '<S330>/Product1'
     */
    rtb_Abs1 = mainV03_56_P.Constant_Value_c - rtb_Rn * rtb_Rn * rtb_Abs1 *
      rtb_Abs1;

    /* Product: '<S328>/Product1' incorporates:
     *  Constant: '<S328>/Constant1'
     *  Sqrt: '<S328>/sqrt'
     */
    rtb_Switch = mainV03_56_P.Constant1_Value_br / sqrt(rtb_Abs1);

    /* Trigonometry: '<S328>/Trigonometric Function1' incorporates:
     *  Constant: '<S328>/Constant'
     *  Constant: '<S328>/Constant2'
     *  Product: '<S328>/Product2'
     *  Product: '<S328>/Product3'
     *  Sum: '<S328>/Sum1'
     */
    mainV03_56_B.TrigonometricFunction1 = rt_atan2d_snf
      (mainV03_56_P.Constant2_Value_l, (mainV03_56_P.Constant_Value_p5 - rtb_Rn *
        rtb_Rn) * rtb_Switch / rtb_Abs1);
  }

  /* Product: '<S313>/x*cos' */
  mainV03_56_B.xcos = rtb_IC4_idx_0 * mainV03_56_B.SinCos_o2;

  /* Product: '<S313>/y*sin' */
  mainV03_56_B.ysin = rtb_IC4_idx_1 * mainV03_56_B.SinCos_o1;

  /* Sum: '<S313>/Sum' */
  mainV03_56_B.Sum = mainV03_56_B.xcos - mainV03_56_B.ysin;

  /* Product: '<S313>/rad lat' */
  mainV03_56_B.radlat = mainV03_56_B.Sum * mainV03_56_B.TrigonometricFunction1;
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Trigonometry: '<S328>/Trigonometric Function2' incorporates:
     *  Constant: '<S328>/Constant3'
     *  Product: '<S328>/Product4'
     *  Trigonometry: '<S328>/Trigonometric Function'
     */
    mainV03_56_B.TrigonometricFunction2 = rt_atan2d_snf
      (mainV03_56_P.Constant3_Value_e, rtb_Switch * cos(rtb_Sum1_no));

    /* Switch: '<S312>/Switch1' incorporates:
     *  Constant: '<S312>/Constant'
     *  Constant: '<S312>/Constant1'
     */
    if (rtb_Compare_e) {
      rtb_Switch = mainV03_56_P.Constant_Value_pw;
    } else {
      rtb_Switch = mainV03_56_P.Constant1_Value_k;
    }

    /* End of Switch: '<S312>/Switch1' */

    /* Sum: '<S312>/Sum' incorporates:
     *  Constant: '<S244>/initial_pos'
     */
    rtb_Sum1_no = rtb_Switch + mainV03_56_P.FlatEarthtoLLA_LL0[1];

    /* Abs: '<S322>/Abs' */
    rtb_Switch = fabs(rtb_Sum1_no);

    /* Switch: '<S322>/Switch' incorporates:
     *  Bias: '<S322>/Bias'
     *  Bias: '<S322>/Bias1'
     *  Constant: '<S322>/Constant2'
     *  Constant: '<S326>/Constant'
     *  Math: '<S322>/Math Function1'
     *  RelationalOperator: '<S326>/Compare'
     */
    if (rtb_Switch > mainV03_56_P.CompareToConstant_const_o) {
      mainV03_56_B.Switch_a = rt_modd_snf(rtb_Sum1_no + mainV03_56_P.Bias_Bias_i,
        mainV03_56_P.Constant2_Value_p) + mainV03_56_P.Bias1_Bias_d;
    } else {
      mainV03_56_B.Switch_a = rtb_Sum1_no;
    }

    /* End of Switch: '<S322>/Switch' */
  }

  /* Product: '<S313>/x*sin' */
  mainV03_56_B.xsin = rtb_IC4_idx_0 * mainV03_56_B.SinCos_o1;

  /* Product: '<S313>/y*cos' */
  mainV03_56_B.ycos = rtb_IC4_idx_1 * mainV03_56_B.SinCos_o2;

  /* Sum: '<S313>/Sum1' */
  mainV03_56_B.Sum1 = mainV03_56_B.xsin + mainV03_56_B.ycos;

  /* Product: '<S313>/rad long ' */
  mainV03_56_B.radlong = mainV03_56_B.TrigonometricFunction2 * mainV03_56_B.Sum1;

  /* Sum: '<S244>/Sum' incorporates:
   *  UnitConversion: '<S314>/Unit Conversion'
   */
  /* Unit Conversion - from: rad to: deg
     Expression: output = (57.2958*input) + (0) */
  mainV03_56_B.Sum_g[0] = 57.295779513082323 * mainV03_56_B.radlat +
    mainV03_56_B.Switch;
  mainV03_56_B.Sum_g[1] = 57.295779513082323 * mainV03_56_B.radlong
    + mainV03_56_B.Switch_a;

  /* Abs: '<S318>/Abs' */
  mainV03_56_B.Abs = fabs(mainV03_56_B.Sum_g[0]);

  /* Switch: '<S318>/Switch' incorporates:
   *  Constant: '<S319>/Constant'
   *  RelationalOperator: '<S319>/Compare'
   */
  if (mainV03_56_B.Abs > mainV03_56_P.CompareToConstant_const_n) {
    /* Bias: '<S318>/Bias' */
    mainV03_56_B.Bias_b = mainV03_56_B.Sum_g[0] + mainV03_56_P.Bias_Bias_g;

    /* Bias: '<S318>/Bias1' incorporates:
     *  Constant: '<S318>/Constant2'
     *  Math: '<S318>/Math Function1'
     */
    mainV03_56_B.Bias1_c = rt_modd_snf(mainV03_56_B.Bias_b,
      mainV03_56_P.Constant2_Value_g) + mainV03_56_P.Bias1_Bias_a;
    mainV03_56_B.Switch_f = mainV03_56_B.Bias1_c;
  } else {
    mainV03_56_B.Switch_f = mainV03_56_B.Sum_g[0];
  }

  /* End of Switch: '<S318>/Switch' */

  /* Abs: '<S315>/Abs1' */
  mainV03_56_B.Abs1 = fabs(mainV03_56_B.Switch_f);

  /* Switch: '<S315>/Switch' incorporates:
   *  Constant: '<S311>/Constant'
   *  Constant: '<S311>/Constant1'
   *  Constant: '<S317>/Constant'
   *  RelationalOperator: '<S317>/Compare'
   *  Switch: '<S311>/Switch1'
   */
  if (mainV03_56_B.Abs1 > mainV03_56_P.CompareToConstant_const_p) {
    /* Bias: '<S315>/Bias' */
    mainV03_56_B.Bias_l = mainV03_56_B.Abs1 + mainV03_56_P.Bias_Bias_n;

    /* Gain: '<S315>/Gain' */
    mainV03_56_B.Gain_o = mainV03_56_P.Gain_Gain_l * mainV03_56_B.Bias_l;

    /* Bias: '<S315>/Bias1' */
    mainV03_56_B.Bias1_g = mainV03_56_B.Gain_o + mainV03_56_P.Bias1_Bias;

    /* Signum: '<S315>/Sign1' */
    if (mainV03_56_B.Switch_f < 0.0) {
      mainV03_56_B.Sign1_f = -1.0;
    } else if (mainV03_56_B.Switch_f > 0.0) {
      mainV03_56_B.Sign1_f = 1.0;
    } else if (mainV03_56_B.Switch_f == 0.0) {
      mainV03_56_B.Sign1_f = 0.0;
    } else {
      mainV03_56_B.Sign1_f = mainV03_56_B.Switch_f;
    }

    /* End of Signum: '<S315>/Sign1' */

    /* Product: '<S315>/Divide1' */
    mainV03_56_B.Divide1_o = mainV03_56_B.Sign1_f * mainV03_56_B.Bias1_g;
    mainV03_56_B.LLA.Latitude_deg = mainV03_56_B.Divide1_o;
    mainV03_56_B.Switch1 = mainV03_56_P.Constant_Value_f;
  } else {
    mainV03_56_B.LLA.Latitude_deg = mainV03_56_B.Switch_f;
    mainV03_56_B.Switch1 = mainV03_56_P.Constant1_Value_a;
  }

  /* End of Switch: '<S315>/Switch' */

  /* Sum: '<S311>/Sum' */
  mainV03_56_B.Sum_c = mainV03_56_B.Switch1 + mainV03_56_B.Sum_g[1];

  /* Abs: '<S316>/Abs' */
  mainV03_56_B.Abs_h = fabs(mainV03_56_B.Sum_c);

  /* Switch: '<S316>/Switch' incorporates:
   *  Constant: '<S320>/Constant'
   *  RelationalOperator: '<S320>/Compare'
   */
  if (mainV03_56_B.Abs_h > mainV03_56_P.CompareToConstant_const_h) {
    /* Bias: '<S316>/Bias' */
    mainV03_56_B.Bias_e = mainV03_56_B.Sum_c + mainV03_56_P.Bias_Bias_e;

    /* Bias: '<S316>/Bias1' incorporates:
     *  Constant: '<S316>/Constant2'
     *  Math: '<S316>/Math Function1'
     */
    mainV03_56_B.Bias1_e = rt_modd_snf(mainV03_56_B.Bias_e,
      mainV03_56_P.Constant2_Value_a) + mainV03_56_P.Bias1_Bias_ae;
    mainV03_56_B.LLA.Longitude_deg = mainV03_56_B.Bias1_e;
  } else {
    mainV03_56_B.LLA.Longitude_deg = mainV03_56_B.Sum_c;
  }

  /* End of Switch: '<S316>/Switch' */

  /* Sum: '<S244>/Sum1' incorporates:
   *  Constant: '<S4>/Constant'
   *  UnaryMinus: '<S244>/Ze2height'
   */
  mainV03_56_B.LLA.Altitude_m = -rtb_IC4_idx_2 - mainV03_56_P.Constant_Value_e;

  /* Integrator: '<S291>/phi theta psi' */
  mainV03_56_B.phithetapsi[0] = mainV03_56_X.phithetapsi_CSTATE[0];
  mainV03_56_B.phithetapsi[1] = mainV03_56_X.phithetapsi_CSTATE[1];
  mainV03_56_B.phithetapsi[2] = mainV03_56_X.phithetapsi_CSTATE[2];

  /* SignalConversion: '<S299>/TmpSignal ConversionAtsincosInport1' */
  mainV03_56_B.TmpSignalConversionAtsincosInport1[0] = mainV03_56_B.phithetapsi
    [2];
  mainV03_56_B.TmpSignalConversionAtsincosInport1[1] = mainV03_56_B.phithetapsi
    [1];
  mainV03_56_B.TmpSignalConversionAtsincosInport1[2] = mainV03_56_B.phithetapsi
    [0];

  /* Trigonometry: '<S299>/sincos' */
  rtb_IC4_idx_0 = sin(mainV03_56_B.TmpSignalConversionAtsincosInport1[0]);
  rtb_IC3_idx_0 = cos(mainV03_56_B.TmpSignalConversionAtsincosInport1[0]);
  rtb_IC4_idx_1 = sin(mainV03_56_B.TmpSignalConversionAtsincosInport1[1]);
  rtb_IC3_idx_1 = cos(mainV03_56_B.TmpSignalConversionAtsincosInport1[1]);
  rtb_IC4_idx_2 = sin(mainV03_56_B.TmpSignalConversionAtsincosInport1[2]);
  rtb_IC3_idx_2 = cos(mainV03_56_B.TmpSignalConversionAtsincosInport1[2]);

  /* Fcn: '<S299>/Fcn11' */
  mainV03_56_B.VectorConcatenate[0] = rtb_IC3_idx_1 * rtb_IC3_idx_0;

  /* Fcn: '<S299>/Fcn21' */
  mainV03_56_B.VectorConcatenate[1] = rtb_IC4_idx_2 * rtb_IC4_idx_1 *
    rtb_IC3_idx_0 - rtb_IC3_idx_2 * rtb_IC4_idx_0;

  /* Fcn: '<S299>/Fcn31' */
  mainV03_56_B.VectorConcatenate[2] = rtb_IC3_idx_2 * rtb_IC4_idx_1 *
    rtb_IC3_idx_0 + rtb_IC4_idx_2 * rtb_IC4_idx_0;

  /* Fcn: '<S299>/Fcn12' */
  mainV03_56_B.VectorConcatenate[3] = rtb_IC3_idx_1 * rtb_IC4_idx_0;

  /* Fcn: '<S299>/Fcn22' */
  mainV03_56_B.VectorConcatenate[4] = rtb_IC4_idx_2 * rtb_IC4_idx_1 *
    rtb_IC4_idx_0 + rtb_IC3_idx_2 * rtb_IC3_idx_0;

  /* Fcn: '<S299>/Fcn32' */
  mainV03_56_B.VectorConcatenate[5] = rtb_IC3_idx_2 * rtb_IC4_idx_1 *
    rtb_IC4_idx_0 - rtb_IC4_idx_2 * rtb_IC3_idx_0;

  /* Fcn: '<S299>/Fcn13' */
  mainV03_56_B.VectorConcatenate[6] = -rtb_IC4_idx_1;

  /* Fcn: '<S299>/Fcn23' */
  mainV03_56_B.VectorConcatenate[7] = rtb_IC4_idx_2 * rtb_IC3_idx_1;

  /* Fcn: '<S299>/Fcn33' */
  mainV03_56_B.VectorConcatenate[8] = rtb_IC3_idx_2 * rtb_IC3_idx_1;
  for (s155_iter = 0; s155_iter < 3; s155_iter++) {
    /* Math: '<S243>/Transpose' */
    mainV03_56_B.Transpose[3 * s155_iter] =
      mainV03_56_B.VectorConcatenate[s155_iter];
    mainV03_56_B.Transpose[1 + 3 * s155_iter] =
      mainV03_56_B.VectorConcatenate[s155_iter + 3];
    mainV03_56_B.Transpose[2 + 3 * s155_iter] =
      mainV03_56_B.VectorConcatenate[s155_iter + 6];

    /* Integrator: '<S243>/ub,vb,wb' */
    mainV03_56_B.ubvbwb[s155_iter] = mainV03_56_X.ubvbwb_CSTATE[s155_iter];
  }

  /* Product: '<S298>/Product' */
  for (rtb_Sum1_eh = 0; rtb_Sum1_eh < 3; rtb_Sum1_eh++) {
    mainV03_56_B.Product[rtb_Sum1_eh] = 0.0;
    mainV03_56_B.Product[rtb_Sum1_eh] += mainV03_56_B.Transpose[rtb_Sum1_eh] *
      mainV03_56_B.ubvbwb[0];
    mainV03_56_B.Product[rtb_Sum1_eh] += mainV03_56_B.Transpose[rtb_Sum1_eh + 3]
      * mainV03_56_B.ubvbwb[1];
    mainV03_56_B.Product[rtb_Sum1_eh] += mainV03_56_B.Transpose[rtb_Sum1_eh + 6]
      * mainV03_56_B.ubvbwb[2];
  }

  /* End of Product: '<S298>/Product' */

  /* InitialCondition: '<S4>/IC11' */
  if ((mainV03_56_DW.IC11_FirstOutputTime == (rtMinusInf)) ||
      (mainV03_56_DW.IC11_FirstOutputTime == mainV03_56_M->Timing.t[0])) {
    mainV03_56_DW.IC11_FirstOutputTime = mainV03_56_M->Timing.t[0];
    rtb_IC3_idx_0 = mainV03_56_P.IC11_Value[0];
    rtb_IC3_idx_1 = mainV03_56_P.IC11_Value[1];
    rtb_IC3_idx_2 = mainV03_56_P.IC11_Value[2];
  } else {
    rtb_IC3_idx_0 = mainV03_56_B.Product[0];
    rtb_IC3_idx_1 = mainV03_56_B.Product[1];
    rtb_IC3_idx_2 = mainV03_56_B.Product[2];
  }

  /* End of InitialCondition: '<S4>/IC11' */

  /* InitialCondition: '<S4>/IC12' */
  if ((mainV03_56_DW.IC12_FirstOutputTime == (rtMinusInf)) ||
      (mainV03_56_DW.IC12_FirstOutputTime == mainV03_56_M->Timing.t[0])) {
    mainV03_56_DW.IC12_FirstOutputTime = mainV03_56_M->Timing.t[0];
    rtb_IC4_idx_0 = mainV03_56_P.IC12_Value[0];
    rtb_IC4_idx_1 = mainV03_56_P.IC12_Value[1];
    rtb_IC4_idx_2 = mainV03_56_P.IC12_Value[2];
  } else {
    rtb_IC4_idx_0 = mainV03_56_B.xeyeze[0];
    rtb_IC4_idx_1 = mainV03_56_B.xeyeze[1];
    rtb_IC4_idx_2 = mainV03_56_B.xeyeze[2];
  }

  /* End of InitialCondition: '<S4>/IC12' */

  /* InitialCondition: '<S4>/IC2' */
  if ((mainV03_56_DW.IC2_FirstOutputTime == (rtMinusInf)) ||
      (mainV03_56_DW.IC2_FirstOutputTime == mainV03_56_M->Timing.t[0])) {
    mainV03_56_DW.IC2_FirstOutputTime = mainV03_56_M->Timing.t[0];
    rtb_sincos_o2_i_idx_0 = mainV03_56_P.IC2_Value_b[0];
    rtb_sincos_o2_i_idx_1 = mainV03_56_P.IC2_Value_b[1];
    rtb_sincos_o2_i_idx_2 = mainV03_56_P.IC2_Value_b[2];
  } else {
    rtb_sincos_o2_i_idx_0 = mainV03_56_B.phithetapsi[0];
    rtb_sincos_o2_i_idx_1 = mainV03_56_B.phithetapsi[1];
    rtb_sincos_o2_i_idx_2 = mainV03_56_B.phithetapsi[2];
  }

  /* End of InitialCondition: '<S4>/IC2' */

  /* BusCreator: '<S4>/EulerBus' */
  mainV03_56_B.Euler.phi = rtb_sincos_o2_i_idx_0;
  mainV03_56_B.Euler.theta = rtb_sincos_o2_i_idx_1;
  mainV03_56_B.Euler.psi = rtb_sincos_o2_i_idx_2;

  /* InitialCondition: '<S4>/IC5' */
  if ((mainV03_56_DW.IC5_FirstOutputTime == (rtMinusInf)) ||
      (mainV03_56_DW.IC5_FirstOutputTime == mainV03_56_M->Timing.t[0])) {
    mainV03_56_DW.IC5_FirstOutputTime = mainV03_56_M->Timing.t[0];
    rtb_sincos_o2_i_idx_0 = mainV03_56_P.IC5_Value[0];
    rtb_sincos_o2_i_idx_1 = mainV03_56_P.IC5_Value[1];
    rtb_sincos_o2_i_idx_2 = mainV03_56_P.IC5_Value[2];
  } else {
    rtb_sincos_o2_i_idx_0 = mainV03_56_B.ubvbwb[0];
    rtb_sincos_o2_i_idx_1 = mainV03_56_B.ubvbwb[1];
    rtb_sincos_o2_i_idx_2 = mainV03_56_B.ubvbwb[2];
  }

  /* End of InitialCondition: '<S4>/IC5' */

  /* Integrator: '<S243>/p,q,r ' */
  mainV03_56_B.pqr[0] = mainV03_56_X.pqr_CSTATE[0];
  mainV03_56_B.pqr[1] = mainV03_56_X.pqr_CSTATE[1];
  mainV03_56_B.pqr[2] = mainV03_56_X.pqr_CSTATE[2];

  /* InitialCondition: '<S4>/IC3' */
  if ((mainV03_56_DW.IC3_FirstOutputTime == (rtMinusInf)) ||
      (mainV03_56_DW.IC3_FirstOutputTime == mainV03_56_M->Timing.t[0])) {
    mainV03_56_DW.IC3_FirstOutputTime = mainV03_56_M->Timing.t[0];
    rtb_sincos_o1_n_idx_0 = mainV03_56_P.IC3_Value_d[0];
    rtb_sincos_o1_n_idx_1 = mainV03_56_P.IC3_Value_d[1];
    rtb_Sum1_no = mainV03_56_P.IC3_Value_d[2];
  } else {
    rtb_sincos_o1_n_idx_0 = mainV03_56_B.pqr[0];
    rtb_sincos_o1_n_idx_1 = mainV03_56_B.pqr[1];
    rtb_Sum1_no = mainV03_56_B.pqr[2];
  }

  /* End of InitialCondition: '<S4>/IC3' */

  /* BusCreator: '<S4>/OmegaBus' */
  mainV03_56_B.Omega_body.p = rtb_sincos_o1_n_idx_0;
  mainV03_56_B.Omega_body.q = rtb_sincos_o1_n_idx_1;
  mainV03_56_B.Omega_body.r = rtb_Sum1_no;

  /* TransferFcn: '<S250>/Transfer Fcn' */
  mainV03_56_B.TransferFcn = 0.0;
  mainV03_56_B.TransferFcn += mainV03_56_P.TransferFcn_C *
    mainV03_56_X.TransferFcn_CSTATE;

  /* TransferFcn: '<S250>/Transfer Fcn1' */
  mainV03_56_B.TransferFcn1 = 0.0;
  mainV03_56_B.TransferFcn1 += mainV03_56_P.TransferFcn1_C *
    mainV03_56_X.TransferFcn1_CSTATE;

  /* TransferFcn: '<S250>/Transfer Fcn2' */
  mainV03_56_B.TransferFcn2 = 0.0;
  mainV03_56_B.TransferFcn2 += mainV03_56_P.TransferFcn2_C *
    mainV03_56_X.TransferFcn2_CSTATE;

  /* SignalConversion: '<S4>/TmpSignal ConversionAtplantDataBusInport8' */
  mainV03_56_B.dOmega_body[0] = mainV03_56_B.TransferFcn;
  mainV03_56_B.dOmega_body[1] = mainV03_56_B.TransferFcn1;
  mainV03_56_B.dOmega_body[2] = mainV03_56_B.TransferFcn2;

  /* TransferFcn: '<S250>/Transfer Fcn3' */
  mainV03_56_B.TransferFcn3 = 0.0;
  mainV03_56_B.TransferFcn3 += mainV03_56_P.TransferFcn3_C *
    mainV03_56_X.TransferFcn3_CSTATE;

  /* TransferFcn: '<S250>/Transfer Fcn4' */
  mainV03_56_B.TransferFcn4 = 0.0;
  mainV03_56_B.TransferFcn4 += mainV03_56_P.TransferFcn4_C *
    mainV03_56_X.TransferFcn4_CSTATE;

  /* TransferFcn: '<S250>/Transfer Fcn5' */
  mainV03_56_B.TransferFcn5 = 0.0;
  mainV03_56_B.TransferFcn5 += mainV03_56_P.TransferFcn5_C *
    mainV03_56_X.TransferFcn5_CSTATE;

  /* SignalConversion: '<S4>/TmpSignal ConversionAtplantDataBusInport9' */
  mainV03_56_B.Accel_body[0] = mainV03_56_B.TransferFcn3;
  mainV03_56_B.Accel_body[1] = mainV03_56_B.TransferFcn4;
  mainV03_56_B.Accel_body[2] = mainV03_56_B.TransferFcn5;

  /* TransferFcn: '<S249>/Transfer Fcn1' */
  mainV03_56_B.TransferFcn1_l = 0.0;
  mainV03_56_B.TransferFcn1_l += mainV03_56_P.TransferFcn1_C_b *
    mainV03_56_X.TransferFcn1_CSTATE_k;

  /* TransferFcn: '<S249>/Transfer Fcn4' */
  mainV03_56_B.TransferFcn4_i = 0.0;
  mainV03_56_B.TransferFcn4_i += mainV03_56_P.TransferFcn4_C_c *
    mainV03_56_X.TransferFcn4_CSTATE_l;

  /* TransferFcn: '<S249>/Transfer Fcn5' */
  mainV03_56_B.TransferFcn5_h = 0.0;
  mainV03_56_B.TransferFcn5_h += mainV03_56_P.TransferFcn5_C_h *
    mainV03_56_X.TransferFcn5_CSTATE_f;

  /* Sum: '<S4>/Sum2' */
  mainV03_56_B.Vaero[0] = rtb_sincos_o2_i_idx_0 - mainV03_56_B.TransferFcn1_l;
  mainV03_56_B.Vaero[1] = rtb_sincos_o2_i_idx_1 - mainV03_56_B.TransferFcn4_i;
  mainV03_56_B.Vaero[2] = rtb_sincos_o2_i_idx_2 - mainV03_56_B.TransferFcn5_h;

  /* InitialCondition: '<S4>/IC7' incorporates:
   *  Trigonometry: '<S245>/Incidence'
   */
  if ((mainV03_56_DW.IC7_FirstOutputTime == (rtMinusInf)) ||
      (mainV03_56_DW.IC7_FirstOutputTime == mainV03_56_M->Timing.t[0])) {
    mainV03_56_DW.IC7_FirstOutputTime = mainV03_56_M->Timing.t[0];
    mainV03_56_B.alpha_j = mainV03_56_P.IC7_Value;
  } else {
    mainV03_56_B.alpha_j = rt_atan2d_snf(mainV03_56_B.Vaero[2],
      mainV03_56_B.Vaero[0]);
  }

  /* End of InitialCondition: '<S4>/IC7' */

  /* Product: '<S335>/Product' */
  mainV03_56_B.Product_c = mainV03_56_B.Vaero[0] * mainV03_56_B.Vaero[0];

  /* Product: '<S335>/Product1' */
  mainV03_56_B.Product1 = mainV03_56_B.Vaero[1] * mainV03_56_B.Vaero[1];

  /* Product: '<S335>/Product2' */
  mainV03_56_B.Product2 = mainV03_56_B.Vaero[2] * mainV03_56_B.Vaero[2];

  /* Sum: '<S335>/Sum' */
  mainV03_56_B.Sum_o = (mainV03_56_B.Product_c + mainV03_56_B.Product1) +
    mainV03_56_B.Product2;

  /* Sqrt: '<S245>/Airspeed' */
  rtb_Sum1_no = sqrt(mainV03_56_B.Sum_o);

  /* Switch: '<S245>/Switch' incorporates:
   *  Constant: '<S245>/Constant'
   */
  if (rtb_Sum1_no > mainV03_56_P.Switch_Threshold_o) {
    /* Product: '<S245>/Product' */
    mainV03_56_B.Product_m = mainV03_56_B.Vaero[1] / rtb_Sum1_no;
    mainV03_56_B.Switch_h = mainV03_56_B.Product_m;
  } else {
    mainV03_56_B.Switch_h = mainV03_56_P.Constant_Value_ep;
  }

  /* End of Switch: '<S245>/Switch' */

  /* InitialCondition: '<S4>/IC8' incorporates:
   *  Trigonometry: '<S245>/Sideslip'
   */
  if ((mainV03_56_DW.IC8_FirstOutputTime == (rtMinusInf)) ||
      (mainV03_56_DW.IC8_FirstOutputTime == mainV03_56_M->Timing.t[0])) {
    mainV03_56_DW.IC8_FirstOutputTime = mainV03_56_M->Timing.t[0];
    mainV03_56_B.beta_l = mainV03_56_P.IC8_Value;
  } else {
    if (mainV03_56_B.Switch_h > 1.0) {
      /* Trigonometry: '<S245>/Sideslip' */
      riseValLimit = 1.0;
    } else if (mainV03_56_B.Switch_h < -1.0) {
      /* Trigonometry: '<S245>/Sideslip' */
      riseValLimit = -1.0;
    } else {
      /* Trigonometry: '<S245>/Sideslip' */
      riseValLimit = mainV03_56_B.Switch_h;
    }

    mainV03_56_B.beta_l = asin(riseValLimit);
  }

  /* End of InitialCondition: '<S4>/IC8' */

  /* Derivative: '<S4>/Derivative' */
  if ((mainV03_56_DW.TimeStampA >= mainV03_56_M->Timing.t[0]) &&
      (mainV03_56_DW.TimeStampB >= mainV03_56_M->Timing.t[0])) {
    mainV03_56_B.Derivative = 0.0;
  } else {
    rtb_Sum1_no = mainV03_56_DW.TimeStampA;
    lastU = &mainV03_56_DW.LastUAtTimeA;
    if (mainV03_56_DW.TimeStampA < mainV03_56_DW.TimeStampB) {
      if (mainV03_56_DW.TimeStampB < mainV03_56_M->Timing.t[0]) {
        rtb_Sum1_no = mainV03_56_DW.TimeStampB;
        lastU = &mainV03_56_DW.LastUAtTimeB;
      }
    } else {
      if (mainV03_56_DW.TimeStampA >= mainV03_56_M->Timing.t[0]) {
        rtb_Sum1_no = mainV03_56_DW.TimeStampB;
        lastU = &mainV03_56_DW.LastUAtTimeB;
      }
    }

    mainV03_56_B.Derivative = (mainV03_56_B.alpha_j - *lastU) /
      (mainV03_56_M->Timing.t[0] - rtb_Sum1_no);
  }

  /* End of Derivative: '<S4>/Derivative' */

  /* RateLimiter: '<S4>/Rate Limiter' */
  if ((mainV03_56_DW.LastMajorTimeA >= mainV03_56_M->Timing.t[0]) &&
      (mainV03_56_DW.LastMajorTimeB >= mainV03_56_M->Timing.t[0])) {
    mainV03_56_B.alpha_dot = mainV03_56_B.Derivative;
  } else {
    if (((mainV03_56_DW.LastMajorTimeA < mainV03_56_DW.LastMajorTimeB) &&
         (mainV03_56_DW.LastMajorTimeB < mainV03_56_M->Timing.t[0])) ||
        ((mainV03_56_DW.LastMajorTimeA >= mainV03_56_DW.LastMajorTimeB) &&
         (mainV03_56_DW.LastMajorTimeA >= mainV03_56_M->Timing.t[0]))) {
      rtb_sincos_o1_n_idx_1 = mainV03_56_M->Timing.t[0] -
        mainV03_56_DW.LastMajorTimeB;
      rtb_Sum1_no = mainV03_56_DW.PrevYB;
    } else {
      rtb_sincos_o1_n_idx_1 = mainV03_56_M->Timing.t[0] -
        mainV03_56_DW.LastMajorTimeA;
      rtb_Sum1_no = mainV03_56_DW.PrevYA;
    }

    riseValLimit = rtb_sincos_o1_n_idx_1 * mainV03_56_P.RateLimiter_RisingLim;
    rtb_sincos_o1_n_idx_0 = mainV03_56_B.Derivative - rtb_Sum1_no;
    if (rtb_sincos_o1_n_idx_0 > riseValLimit) {
      mainV03_56_B.alpha_dot = rtb_Sum1_no + riseValLimit;
    } else {
      rtb_sincos_o1_n_idx_1 *= mainV03_56_P.RateLimiter_FallingLim;
      if (rtb_sincos_o1_n_idx_0 < rtb_sincos_o1_n_idx_1) {
        mainV03_56_B.alpha_dot = rtb_Sum1_no + rtb_sincos_o1_n_idx_1;
      } else {
        mainV03_56_B.alpha_dot = mainV03_56_B.Derivative;
      }
    }
  }

  /* End of RateLimiter: '<S4>/Rate Limiter' */

  /* Derivative: '<S4>/Derivative1' */
  if ((mainV03_56_DW.TimeStampA_i >= mainV03_56_M->Timing.t[0]) &&
      (mainV03_56_DW.TimeStampB_a >= mainV03_56_M->Timing.t[0])) {
    mainV03_56_B.Derivative1 = 0.0;
  } else {
    rtb_Sum1_no = mainV03_56_DW.TimeStampA_i;
    lastU = &mainV03_56_DW.LastUAtTimeA_o;
    if (mainV03_56_DW.TimeStampA_i < mainV03_56_DW.TimeStampB_a) {
      if (mainV03_56_DW.TimeStampB_a < mainV03_56_M->Timing.t[0]) {
        rtb_Sum1_no = mainV03_56_DW.TimeStampB_a;
        lastU = &mainV03_56_DW.LastUAtTimeB_j;
      }
    } else {
      if (mainV03_56_DW.TimeStampA_i >= mainV03_56_M->Timing.t[0]) {
        rtb_Sum1_no = mainV03_56_DW.TimeStampB_a;
        lastU = &mainV03_56_DW.LastUAtTimeB_j;
      }
    }

    mainV03_56_B.Derivative1 = (mainV03_56_B.beta_l - *lastU) /
      (mainV03_56_M->Timing.t[0] - rtb_Sum1_no);
  }

  /* End of Derivative: '<S4>/Derivative1' */

  /* RateLimiter: '<S4>/Rate Limiter1' */
  if ((mainV03_56_DW.LastMajorTimeA_f >= mainV03_56_M->Timing.t[0]) &&
      (mainV03_56_DW.LastMajorTimeB_d >= mainV03_56_M->Timing.t[0])) {
    mainV03_56_B.beta_dot = mainV03_56_B.Derivative1;
  } else {
    if (((mainV03_56_DW.LastMajorTimeA_f < mainV03_56_DW.LastMajorTimeB_d) &&
         (mainV03_56_DW.LastMajorTimeB_d < mainV03_56_M->Timing.t[0])) ||
        ((mainV03_56_DW.LastMajorTimeA_f >= mainV03_56_DW.LastMajorTimeB_d) &&
         (mainV03_56_DW.LastMajorTimeA_f >= mainV03_56_M->Timing.t[0]))) {
      rtb_sincos_o1_n_idx_1 = mainV03_56_M->Timing.t[0] -
        mainV03_56_DW.LastMajorTimeB_d;
      rtb_Sum1_no = mainV03_56_DW.PrevYB_d;
    } else {
      rtb_sincos_o1_n_idx_1 = mainV03_56_M->Timing.t[0] -
        mainV03_56_DW.LastMajorTimeA_f;
      rtb_Sum1_no = mainV03_56_DW.PrevYA_o;
    }

    riseValLimit = rtb_sincos_o1_n_idx_1 * mainV03_56_P.RateLimiter1_RisingLim;
    rtb_sincos_o1_n_idx_0 = mainV03_56_B.Derivative1 - rtb_Sum1_no;
    if (rtb_sincos_o1_n_idx_0 > riseValLimit) {
      mainV03_56_B.beta_dot = rtb_Sum1_no + riseValLimit;
    } else {
      rtb_sincos_o1_n_idx_1 *= mainV03_56_P.RateLimiter1_FallingLim;
      if (rtb_sincos_o1_n_idx_0 < rtb_sincos_o1_n_idx_1) {
        mainV03_56_B.beta_dot = rtb_Sum1_no + rtb_sincos_o1_n_idx_1;
      } else {
        mainV03_56_B.beta_dot = mainV03_56_B.Derivative1;
      }
    }
  }

  /* End of RateLimiter: '<S4>/Rate Limiter1' */
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* SignalConversion: '<S336>/TmpSignal ConversionAtPreLook-Up Index SearchInport2' incorporates:
     *  Constant: '<S336>/Constant'
     *  Constant: '<S336>/Constant1'
     */
    rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[0] =
      mainV03_56_P.EstimateCenterofGravity_emass;
    rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[1] =
      mainV03_56_P.EstimateCenterofGravity_fmass;

    /* PreLookup: '<S336>/PreLook-Up Index Search' incorporates:
     *  Constant: '<S246>/mass (kg)'
     */
    rtb_PreLookUpIndexSearch_o1 = plook_s32dd_bincp(mainV03_56_P.masskg_Value,
      rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2, 1U, &rtb_Switch,
      &mainV03_56_DW.PreLookUpIndexSearch_DWORK1);

    /* LookupNDDirect: '<S340>/[0]'
     *
     * About '<S340>/[0]':
     *  3-dimensional Direct Look-Up returning a 2-D Matrix
     */
    if (rtb_PreLookUpIndexSearch_o1 > 1) {
      s106_iter = 1;
    } else if (rtb_PreLookUpIndexSearch_o1 < 0) {
      s106_iter = 0;
    } else {
      s106_iter = rtb_PreLookUpIndexSearch_o1;
    }

    s106_iter *= 3;

    /* Sum: '<S340>/Sum' incorporates:
     *  Constant: '<S340>/offset for upper index 1'
     */
    if ((rtb_PreLookUpIndexSearch_o1 < 0) &&
        (mainV03_56_P.offsetforupperindex1_Value < MIN_int32_T
         - rtb_PreLookUpIndexSearch_o1)) {
      rtb_PreLookUpIndexSearch_o1 = MIN_int32_T;
    } else if ((rtb_PreLookUpIndexSearch_o1 > 0) &&
               (mainV03_56_P.offsetforupperindex1_Value > MAX_int32_T
                - rtb_PreLookUpIndexSearch_o1)) {
      rtb_PreLookUpIndexSearch_o1 = MAX_int32_T;
    } else {
      rtb_PreLookUpIndexSearch_o1 += mainV03_56_P.offsetforupperindex1_Value;
    }

    /* LookupNDDirect: '<S340>/[1]' incorporates:
     *  Sum: '<S340>/Sum'
     *
     * About '<S340>/[1]':
     *  3-dimensional Direct Look-Up returning a 2-D Matrix
     */
    if (rtb_PreLookUpIndexSearch_o1 > 1) {
      rtb_PreLookUpIndexSearch_o1 = 1;
    } else {
      if (rtb_PreLookUpIndexSearch_o1 < 0) {
        rtb_PreLookUpIndexSearch_o1 = 0;
      }
    }

    rtb_PreLookUpIndexSearch_o1 *= 3;

    /* Sum: '<S340>/Sum3' incorporates:
     *  Fcn: '<S340>/Fcn'
     *  LookupNDDirect: '<S340>/[0]'
     *  LookupNDDirect: '<S340>/[1]'
     *  Product: '<S340>/1-lambda_x'
     *  Product: '<S340>/lambda_x.'
     *
     * About '<S340>/[0]':
     *  3-dimensional Direct Look-Up returning a 2-D Matrix
     *
     * About '<S340>/[1]':
     *  3-dimensional Direct Look-Up returning a 2-D Matrix
     */
    mainV03_56_B.Sum3[0] = (1.0 - rtb_Switch) *
      mainV03_56_P.InterpolateCG_matrix[s106_iter] +
      mainV03_56_P.InterpolateCG_matrix[rtb_PreLookUpIndexSearch_o1] *
      rtb_Switch;
    mainV03_56_B.Sum3[1] = mainV03_56_P.InterpolateCG_matrix[1 + s106_iter] *
      (1.0 - rtb_Switch) + mainV03_56_P.InterpolateCG_matrix[1 +
      rtb_PreLookUpIndexSearch_o1] * rtb_Switch;
    mainV03_56_B.Sum3[2] = mainV03_56_P.InterpolateCG_matrix[2 + s106_iter] *
      (1.0 - rtb_Switch) + mainV03_56_P.InterpolateCG_matrix[2 +
      rtb_PreLookUpIndexSearch_o1] * rtb_Switch;

    /* Product: '<S247>/Product2' incorporates:
     *  Constant: '<S343>/Battery_Capacity'
     *  Constant: '<S343>/profundidad_descarga'
     */
    mainV03_56_B.nominalCapacityAh = mainV03_56_P.Battery_Capacity_Value *
      mainV03_56_P.profundidad_descarga_Value;
  }

  /* Integrator: '<S247>/Integrator' */
  mainV03_56_B.usedCapacityAs = mainV03_56_X.Integrator_CSTATE;

  /* Product: '<S247>/Divide' incorporates:
   *  Constant: '<S247>/Constant'
   */
  mainV03_56_B.usedCapacityAh = mainV03_56_B.usedCapacityAs /
    mainV03_56_P.Constant_Value_co;

  /* Sum: '<S247>/Sum6' */
  mainV03_56_B.remainingCapacityAh = mainV03_56_B.nominalCapacityAh -
    mainV03_56_B.usedCapacityAh;

  /* BusCreator: '<S4>/plantDataBus' */
  mainV03_56_B.plantData.LLA = mainV03_56_B.LLA;
  mainV03_56_B.plantData.V_ned[0] = rtb_IC3_idx_0;
  mainV03_56_B.plantData.X_ned[0] = rtb_IC4_idx_0;
  mainV03_56_B.plantData.V_ned[1] = rtb_IC3_idx_1;
  mainV03_56_B.plantData.X_ned[1] = rtb_IC4_idx_1;
  mainV03_56_B.plantData.V_ned[2] = rtb_IC3_idx_2;
  mainV03_56_B.plantData.X_ned[2] = rtb_IC4_idx_2;
  mainV03_56_B.plantData.Euler = mainV03_56_B.Euler;
  memcpy(&mainV03_56_B.plantData.DCM_body_earth[0],
         &mainV03_56_B.VectorConcatenate[0], 9U * sizeof(real_T));
  mainV03_56_B.plantData.Omega_body = mainV03_56_B.Omega_body;
  mainV03_56_B.plantData.alpha = mainV03_56_B.alpha_j;
  mainV03_56_B.plantData.beta = mainV03_56_B.beta_l;
  mainV03_56_B.plantData.alpha_dot = mainV03_56_B.alpha_dot;
  mainV03_56_B.plantData.beta_dot = mainV03_56_B.beta_dot;
  mainV03_56_B.plantData.V_body[0] = rtb_sincos_o2_i_idx_0;
  mainV03_56_B.plantData.dOmega_body[0] = mainV03_56_B.dOmega_body[0];
  mainV03_56_B.plantData.Accel_body[0] = mainV03_56_B.Accel_body[0];
  mainV03_56_B.plantData.CG[0] = mainV03_56_B.Sum3[0];
  mainV03_56_B.plantData.V_body[1] = rtb_sincos_o2_i_idx_1;
  mainV03_56_B.plantData.dOmega_body[1] = mainV03_56_B.dOmega_body[1];
  mainV03_56_B.plantData.Accel_body[1] = mainV03_56_B.Accel_body[1];
  mainV03_56_B.plantData.CG[1] = mainV03_56_B.Sum3[1];
  mainV03_56_B.plantData.V_body[2] = rtb_sincos_o2_i_idx_2;
  mainV03_56_B.plantData.dOmega_body[2] = mainV03_56_B.dOmega_body[2];
  mainV03_56_B.plantData.Accel_body[2] = mainV03_56_B.Accel_body[2];
  mainV03_56_B.plantData.CG[2] = mainV03_56_B.Sum3[2];
  mainV03_56_B.plantData.remainingCapacity = mainV03_56_B.remainingCapacityAh;

  /* SignalConversion: '<S7>/SigConversion_InsertedFor_Bus Selector4_at_outport_0' */
  mainV03_56_B.alpha_f = mainV03_56_B.plantData.alpha;

  /* SignalConversion: '<S7>/SigConversion_InsertedFor_Bus Selector4_at_outport_1' */
  mainV03_56_B.beta_n = mainV03_56_B.plantData.beta;
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
  }

  /* SignalConversion: '<S7>/SigConversion_InsertedFor_Bus Selector2_at_outport_0' */
  mainV03_56_B.Altitude_m = mainV03_56_B.plantData.LLA.Altitude_m;
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
  }

  /* SignalConversion: '<S7>/SigConversion_InsertedFor_Bus Selector_at_outport_0' */
  mainV03_56_B.X_ned[0] = mainV03_56_B.plantData.X_ned[0];
  mainV03_56_B.X_ned[1] = mainV03_56_B.plantData.X_ned[1];
  mainV03_56_B.X_ned[2] = mainV03_56_B.plantData.X_ned[2];
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
  }

  /* SignalConversion: '<S7>/SigConversion_InsertedFor_Bus Selector1_at_outport_0' */
  mainV03_56_B.phi = mainV03_56_B.plantData.Euler.phi;

  /* SignalConversion: '<S7>/SigConversion_InsertedFor_Bus Selector1_at_outport_1' */
  mainV03_56_B.theta = mainV03_56_B.plantData.Euler.theta;

  /* SignalConversion: '<S7>/SigConversion_InsertedFor_Bus Selector1_at_outport_2' */
  mainV03_56_B.psi = mainV03_56_B.plantData.Euler.psi;
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
  }

  /* SignalConversion: '<S7>/SigConversion_InsertedFor_Bus Selector5_at_outport_0' */
  mainV03_56_B.Latitude_deg = mainV03_56_B.plantData.LLA.Latitude_deg;

  /* SignalConversion: '<S7>/SigConversion_InsertedFor_Bus Selector5_at_outport_1' */
  mainV03_56_B.Longitude_deg = mainV03_56_B.plantData.LLA.Longitude_deg;
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
  }

  /* SignalConversion: '<S552>/TmpSignal ConversionAtsincosInport1' incorporates:
   *  Gain: '<S7>/Gain'
   */
  rtb_sincos_o1_n_idx_0 = mainV03_56_P.Gain_Gain_e * mainV03_56_B.phi;

  /* Trigonometry: '<S552>/sincos' incorporates:
   *  SignalConversion: '<S552>/TmpSignal ConversionAtsincosInport1'
   */
  rtb_sincos_o2_i_idx_0 = cos(rtb_sincos_o1_n_idx_0);
  rtb_sincos_o1_n_idx_0 = sin(rtb_sincos_o1_n_idx_0);
  rtb_sincos_o2_i_idx_1 = cos(mainV03_56_B.theta);
  rtb_sincos_o1_n_idx_1 = sin(mainV03_56_B.theta);
  rtb_Sum1_no = cos(mainV03_56_B.psi);
  rtb_IC4_idx_0 = sin(mainV03_56_B.psi);

  /* Fcn: '<S552>/Fcn11' incorporates:
   *  Trigonometry: '<S552>/sincos'
   */
  rtb_IC3_idx_1 = rtb_Sum1_no * rtb_sincos_o2_i_idx_0 - rtb_sincos_o1_n_idx_1 *
    rtb_IC4_idx_0 * rtb_sincos_o1_n_idx_0;

  /* Fcn: '<S552>/Fcn12' incorporates:
   *  Trigonometry: '<S552>/sincos'
   */
  rtb_IC4_idx_1 = rtb_sincos_o1_n_idx_1 * rtb_IC4_idx_0 * rtb_sincos_o2_i_idx_0
    + rtb_Sum1_no * rtb_sincos_o1_n_idx_0;

  /* Fcn: '<S552>/Fcn22' */
  rtb_IC3_idx_2 = rtb_sincos_o2_i_idx_1 * rtb_sincos_o2_i_idx_0;

  /* Fcn: '<S552>/Fcn13' */
  rtb_IC3_idx_0 = -rtb_IC4_idx_0 * rtb_sincos_o2_i_idx_1;

  /* Fcn: '<S552>/Fcn33' incorporates:
   *  Trigonometry: '<S552>/sincos'
   */
  rtb_sincos_o2_i_idx_2 = rtb_sincos_o2_i_idx_1 * rtb_Sum1_no;

  /* Sum: '<S553>/Sum of Elements3' */
  rtb_IC4_idx_2 = (rtb_IC3_idx_1 + rtb_IC3_idx_2) + rtb_sincos_o2_i_idx_2;

  /* If: '<S553>/If' incorporates:
   *  Abs: '<S553>/Abs'
   *  Constant: '<S553>/Constant'
   *  Sum: '<S553>/Sum'
   */
  if (rtmIsMajorTimeStep(mainV03_56_M)) {
    if (rtb_IC4_idx_2 + mainV03_56_P.Constant_Value_fu[0] >= 3.0) {
      rtAction = 0;
    } else if (fabs(rtb_IC4_idx_2 + mainV03_56_P.Constant_Value_fu[1]) <=
               1.0E-12) {
      rtAction = 1;
    } else {
      rtAction = 2;
    }

    mainV03_56_DW.If_ActiveSubsystem = rtAction;
  } else {
    rtAction = mainV03_56_DW.If_ActiveSubsystem;
  }

  switch (rtAction) {
   case 0:
    /* Outputs for IfAction SubSystem: '<S553>/Phi == 0' incorporates:
     *  ActionPort: '<S557>/Action Port'
     */
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* SignalConversion: '<S557>/OutportBufferForVRRot_0' incorporates:
       *  Constant: '<S557>/Trace=3=>Phi=0'
       */
      mainV03_56_B.Merge[0] = mainV03_56_P.Trace3Phi0_Value[0];
      mainV03_56_B.Merge[1] = mainV03_56_P.Trace3Phi0_Value[1];
      mainV03_56_B.Merge[2] = mainV03_56_P.Trace3Phi0_Value[2];
      mainV03_56_B.Merge[3] = mainV03_56_P.Trace3Phi0_Value[3];
    }

    /* End of Outputs for SubSystem: '<S553>/Phi == 0' */
    break;

   case 1:
    /* Outputs for IfAction SubSystem: '<S553>/Phi == pi' incorporates:
     *  ActionPort: '<S558>/Action Port'
     */
    /* Sum: '<S558>/Sum' incorporates:
     *  Constant: '<S558>/Constant'
     */
    mainV03_56_B.Sum_bt[0] = mainV03_56_P.Constant_Value_hj + rtb_IC3_idx_1;
    mainV03_56_B.Sum_bt[1] = mainV03_56_P.Constant_Value_hj + rtb_IC3_idx_2;
    mainV03_56_B.Sum_bt[2] = mainV03_56_P.Constant_Value_hj +
      rtb_sincos_o2_i_idx_2;

    /* DeadZone: '<S559>/Dead Zone' incorporates:
     *  Fcn: '<S552>/Fcn23'
     */
    if (rtb_sincos_o1_n_idx_1 >
        mainV03_56_P.RotationMatrixtoVRMLRotation_maxzero) {
      mainV03_56_B.DeadZone[0] = rtb_sincos_o1_n_idx_1 -
        mainV03_56_P.RotationMatrixtoVRMLRotation_maxzero;
    } else if (rtb_sincos_o1_n_idx_1 >=
               -mainV03_56_P.RotationMatrixtoVRMLRotation_maxzero) {
      mainV03_56_B.DeadZone[0] = 0.0;
    } else {
      mainV03_56_B.DeadZone[0] = rtb_sincos_o1_n_idx_1 -
        (-mainV03_56_P.RotationMatrixtoVRMLRotation_maxzero);
    }

    if (rtb_IC3_idx_0 > mainV03_56_P.RotationMatrixtoVRMLRotation_maxzero) {
      mainV03_56_B.DeadZone[1] = rtb_IC3_idx_0 -
        mainV03_56_P.RotationMatrixtoVRMLRotation_maxzero;
    } else if (rtb_IC3_idx_0 >=
               -mainV03_56_P.RotationMatrixtoVRMLRotation_maxzero) {
      mainV03_56_B.DeadZone[1] = 0.0;
    } else {
      mainV03_56_B.DeadZone[1] = rtb_IC3_idx_0 -
        (-mainV03_56_P.RotationMatrixtoVRMLRotation_maxzero);
    }

    if (rtb_IC4_idx_1 > mainV03_56_P.RotationMatrixtoVRMLRotation_maxzero) {
      mainV03_56_B.DeadZone[2] = rtb_IC4_idx_1 -
        mainV03_56_P.RotationMatrixtoVRMLRotation_maxzero;
    } else if (rtb_IC4_idx_1 >=
               -mainV03_56_P.RotationMatrixtoVRMLRotation_maxzero) {
      mainV03_56_B.DeadZone[2] = 0.0;
    } else {
      mainV03_56_B.DeadZone[2] = rtb_IC4_idx_1 -
        (-mainV03_56_P.RotationMatrixtoVRMLRotation_maxzero);
    }

    /* End of DeadZone: '<S559>/Dead Zone' */

    /* Gain: '<S558>/Gain' */
    mainV03_56_B.Gain_c[0] = mainV03_56_P.Gain_Gain_a * mainV03_56_B.Sum_bt[0];

    /* Signum: '<S559>/Sign1' */
    if (mainV03_56_B.DeadZone[0] < 0.0) {
      mainV03_56_B.Sign1[0] = -1.0;
    } else if (mainV03_56_B.DeadZone[0] > 0.0) {
      mainV03_56_B.Sign1[0] = 1.0;
    } else if (mainV03_56_B.DeadZone[0] == 0.0) {
      mainV03_56_B.Sign1[0] = 0.0;
    } else {
      mainV03_56_B.Sign1[0] = mainV03_56_B.DeadZone[0];
    }

    /* Gain: '<S558>/Gain' */
    mainV03_56_B.Gain_c[1] = mainV03_56_P.Gain_Gain_a * mainV03_56_B.Sum_bt[1];

    /* Signum: '<S559>/Sign1' */
    if (mainV03_56_B.DeadZone[1] < 0.0) {
      mainV03_56_B.Sign1[1] = -1.0;
    } else if (mainV03_56_B.DeadZone[1] > 0.0) {
      mainV03_56_B.Sign1[1] = 1.0;
    } else if (mainV03_56_B.DeadZone[1] == 0.0) {
      mainV03_56_B.Sign1[1] = 0.0;
    } else {
      mainV03_56_B.Sign1[1] = mainV03_56_B.DeadZone[1];
    }

    /* Gain: '<S558>/Gain' */
    mainV03_56_B.Gain_c[2] = mainV03_56_P.Gain_Gain_a * mainV03_56_B.Sum_bt[2];

    /* Signum: '<S559>/Sign1' */
    if (mainV03_56_B.DeadZone[2] < 0.0) {
      mainV03_56_B.Sign1[2] = -1.0;
    } else if (mainV03_56_B.DeadZone[2] > 0.0) {
      mainV03_56_B.Sign1[2] = 1.0;
    } else if (mainV03_56_B.DeadZone[2] == 0.0) {
      mainV03_56_B.Sign1[2] = 0.0;
    } else {
      mainV03_56_B.Sign1[2] = mainV03_56_B.DeadZone[2];
    }

    /* Sum: '<S559>/Sum of Elements3' */
    mainV03_56_B.SumofElements3 = (mainV03_56_B.Sign1[0] + mainV03_56_B.Sign1[1])
      + mainV03_56_B.Sign1[2];

    /* Switch: '<S559>/Switch' incorporates:
     *  Constant: '<S559>/Pi1'
     */
    if (mainV03_56_B.SumofElements3 >= mainV03_56_P.Switch_Threshold_p) {
      mainV03_56_B.Switch_dl[0] = mainV03_56_P.Pi1_Value[0];
      mainV03_56_B.Switch_dl[1] = mainV03_56_P.Pi1_Value[1];
      mainV03_56_B.Switch_dl[2] = mainV03_56_P.Pi1_Value[2];
    } else {
      /* Product: '<S559>/Product of Elements' */
      mainV03_56_B.ProductofElements = mainV03_56_B.Sign1[0] *
        mainV03_56_B.Sign1[1] * mainV03_56_B.Sign1[2];

      /* Switch: '<S559>/Switch1' incorporates:
       *  CombinatorialLogic: '<S559>/Shift right'
       */
      if (mainV03_56_B.ProductofElements != 0.0) {
        /* Gain: '<S559>/Gain1' */
        mainV03_56_B.Gain1_m[0] = mainV03_56_P.Gain1_Gain_kf *
          mainV03_56_B.Sign1[0];
        mainV03_56_B.Switch1_a[0] = mainV03_56_B.Gain1_m[0];

        /* Gain: '<S559>/Gain1' */
        mainV03_56_B.Gain1_m[1] = mainV03_56_P.Gain1_Gain_kf *
          mainV03_56_B.Sign1[1];
        mainV03_56_B.Switch1_a[1] = mainV03_56_B.Gain1_m[1];

        /* Gain: '<S559>/Gain1' */
        mainV03_56_B.Gain1_m[2] = mainV03_56_P.Gain1_Gain_kf *
          mainV03_56_B.Sign1[2];
        mainV03_56_B.Switch1_a[2] = mainV03_56_B.Gain1_m[2];
      } else {
        /* CombinatorialLogic: '<S559>/Shift right' incorporates:
         *  Constant: '<S560>/Constant'
         *  RelationalOperator: '<S560>/Compare'
         */
        s106_iter = (int32_T)(((((uint32_T)(mainV03_56_B.Sign1[0] !=
          mainV03_56_P.Constant_Value_jg) << 1) + (mainV03_56_B.Sign1[1] !=
          mainV03_56_P.Constant_Value_jg)) << 1) + (mainV03_56_B.Sign1[2] !=
          mainV03_56_P.Constant_Value_jg));
        mainV03_56_B.Switch1_a[0] = mainV03_56_P.Shiftright_table[(uint32_T)
          s106_iter];
        mainV03_56_B.Switch1_a[1] = mainV03_56_P.Shiftright_table[s106_iter + 8U];
        mainV03_56_B.Switch1_a[2] = mainV03_56_P.Shiftright_table[s106_iter +
          16U];
      }

      /* End of Switch: '<S559>/Switch1' */
      mainV03_56_B.Switch_dl[0] = mainV03_56_B.Switch1_a[0];
      mainV03_56_B.Switch_dl[1] = mainV03_56_B.Switch1_a[1];
      mainV03_56_B.Switch_dl[2] = mainV03_56_B.Switch1_a[2];
    }

    /* End of Switch: '<S559>/Switch' */

    /* MinMax: '<S558>/MinMax' */
    if (mainV03_56_B.Gain_c[0] >= 0.0) {
      mainV03_56_B.MinMax[0] = mainV03_56_B.Gain_c[0];
    } else {
      mainV03_56_B.MinMax[0] = 0.0;
    }

    /* Math: '<S558>/Math Function'
     *
     * About '<S558>/Math Function':
     *  Operator: sqrt
     */
    if (mainV03_56_B.MinMax[0] < 0.0) {
      rtb_Sum1_no = -sqrt(fabs(mainV03_56_B.MinMax[0]));
    } else {
      rtb_Sum1_no = sqrt(mainV03_56_B.MinMax[0]);
    }

    /* Switch: '<S558>/Switch' */
    if (rtb_Sum1_no > mainV03_56_P.RotationMatrixtoVRMLRotation_maxzero) {
      mainV03_56_B.Switch_dj[0] = rtb_Sum1_no;
    } else {
      mainV03_56_B.Switch_dj[0] = 0.0;
    }

    /* Product: '<S558>/Product' */
    mainV03_56_B.Product_o[0] = mainV03_56_B.Switch_dj[0] *
      mainV03_56_B.Switch_dl[0];

    /* SignalConversion: '<S558>/OutportBufferForVRRot_Pi' */
    mainV03_56_B.Merge[0] = mainV03_56_B.Product_o[0];

    /* MinMax: '<S558>/MinMax' */
    if (mainV03_56_B.Gain_c[1] >= 0.0) {
      mainV03_56_B.MinMax[1] = mainV03_56_B.Gain_c[1];
    } else {
      mainV03_56_B.MinMax[1] = 0.0;
    }

    /* Math: '<S558>/Math Function'
     *
     * About '<S558>/Math Function':
     *  Operator: sqrt
     */
    if (mainV03_56_B.MinMax[1] < 0.0) {
      rtb_Sum1_no = -sqrt(fabs(mainV03_56_B.MinMax[1]));
    } else {
      rtb_Sum1_no = sqrt(mainV03_56_B.MinMax[1]);
    }

    /* Switch: '<S558>/Switch' */
    if (rtb_Sum1_no > mainV03_56_P.RotationMatrixtoVRMLRotation_maxzero) {
      mainV03_56_B.Switch_dj[1] = rtb_Sum1_no;
    } else {
      mainV03_56_B.Switch_dj[1] = 0.0;
    }

    /* Product: '<S558>/Product' */
    mainV03_56_B.Product_o[1] = mainV03_56_B.Switch_dj[1] *
      mainV03_56_B.Switch_dl[1];

    /* SignalConversion: '<S558>/OutportBufferForVRRot_Pi' */
    mainV03_56_B.Merge[1] = mainV03_56_B.Product_o[1];

    /* MinMax: '<S558>/MinMax' */
    if (mainV03_56_B.Gain_c[2] >= 0.0) {
      mainV03_56_B.MinMax[2] = mainV03_56_B.Gain_c[2];
    } else {
      mainV03_56_B.MinMax[2] = 0.0;
    }

    /* Math: '<S558>/Math Function'
     *
     * About '<S558>/Math Function':
     *  Operator: sqrt
     */
    if (mainV03_56_B.MinMax[2] < 0.0) {
      rtb_Sum1_no = -sqrt(fabs(mainV03_56_B.MinMax[2]));
    } else {
      rtb_Sum1_no = sqrt(mainV03_56_B.MinMax[2]);
    }

    /* Switch: '<S558>/Switch' */
    if (rtb_Sum1_no > mainV03_56_P.RotationMatrixtoVRMLRotation_maxzero) {
      mainV03_56_B.Switch_dj[2] = rtb_Sum1_no;
    } else {
      mainV03_56_B.Switch_dj[2] = 0.0;
    }

    /* Product: '<S558>/Product' */
    mainV03_56_B.Product_o[2] = mainV03_56_B.Switch_dj[2] *
      mainV03_56_B.Switch_dl[2];

    /* SignalConversion: '<S558>/OutportBufferForVRRot_Pi' incorporates:
     *  Constant: '<S558>/Pi'
     */
    mainV03_56_B.Merge[2] = mainV03_56_B.Product_o[2];
    mainV03_56_B.Merge[3] = mainV03_56_P.Pi_Value;

    /* End of Outputs for SubSystem: '<S553>/Phi == pi' */
    break;

   case 2:
    /* Outputs for IfAction SubSystem: '<S553>/General case' incorporates:
     *  ActionPort: '<S556>/Action Port'
     */
    /* Sum: '<S556>/Subtract' incorporates:
     *  Constant: '<S556>/Constant'
     */
    mainV03_56_B.Subtract = rtb_IC4_idx_2 - mainV03_56_P.Constant_Value_nm;

    /* Gain: '<S556>/Gain1' */
    mainV03_56_B.Gain1 = mainV03_56_P.Gain1_Gain_a * mainV03_56_B.Subtract;

    /* Trigonometry: '<S556>/Trigonometric Function1' */
    if (mainV03_56_B.Gain1 > 1.0) {
      riseValLimit = 1.0;
    } else if (mainV03_56_B.Gain1 < -1.0) {
      riseValLimit = -1.0;
    } else {
      riseValLimit = mainV03_56_B.Gain1;
    }

    rtb_IC4_idx_2 = acos(riseValLimit);

    /* End of Trigonometry: '<S556>/Trigonometric Function1' */

    /* DataTypeConversion: '<S556>/Data Type Conversion' */
    mainV03_56_B.Merge[3] = rtb_IC4_idx_2;

    /* Sum: '<S556>/Sum' incorporates:
     *  Fcn: '<S552>/Fcn21'
     *  Fcn: '<S552>/Fcn23'
     *  Fcn: '<S552>/Fcn31'
     *  Fcn: '<S552>/Fcn32'
     *  Trigonometry: '<S552>/sincos'
     */
    mainV03_56_B.Sum_eh[0] = (rtb_IC4_idx_0 * rtb_sincos_o1_n_idx_0 -
      rtb_sincos_o1_n_idx_1 * rtb_Sum1_no * rtb_sincos_o2_i_idx_0) -
      rtb_sincos_o1_n_idx_1;
    mainV03_56_B.Sum_eh[1] = rtb_IC3_idx_0 - (rtb_sincos_o1_n_idx_1 *
      rtb_Sum1_no * rtb_sincos_o1_n_idx_0 + rtb_IC4_idx_0 *
      rtb_sincos_o2_i_idx_0);
    mainV03_56_B.Sum_eh[2] = -rtb_sincos_o2_i_idx_1 * rtb_sincos_o1_n_idx_0 -
      rtb_IC4_idx_1;

    /* Gain: '<S556>/Gain' incorporates:
     *  Trigonometry: '<S556>/Trigonometric Function'
     */
    mainV03_56_B.Gain_a = mainV03_56_P.Gain_Gain_ao * sin(rtb_IC4_idx_2);

    /* Product: '<S556>/Divide' */
    mainV03_56_B.Merge[0] = mainV03_56_B.Sum_eh[0] / mainV03_56_B.Gain_a;
    mainV03_56_B.Merge[1] = mainV03_56_B.Sum_eh[1] / mainV03_56_B.Gain_a;
    mainV03_56_B.Merge[2] = mainV03_56_B.Sum_eh[2] / mainV03_56_B.Gain_a;

    /* End of Outputs for SubSystem: '<S553>/General case' */
    break;
  }

  /* End of If: '<S553>/If' */

  /* Gain: '<S551>/rad-->deg' */
  rtb_sincos_o2_p[0] = mainV03_56_P.raddeg_Gain * mainV03_56_B.X_ned[1];

  /* Gain: '<S551>/rad-->deg1' incorporates:
   *  Constant: '<S551>/Constant'
   *  Sum: '<S551>/Sum'
   */
  rtb_sincos_o2_p[1] = (mainV03_56_B.X_ned[2] + mainV03_56_P.Constant_Value_jm) *
    mainV03_56_P.raddeg1_Gain;

  /* SignalConversion: '<S551>/ConcatBufferAtConcatenateIn3' */
  rtb_sincos_o2_p[2] = mainV03_56_B.X_ned[0];

  /* SignalConversion: '<S7>/SigConversion_InsertedFor_Bus Selector3_at_outport_0' */
  mainV03_56_B.V_body[0] = mainV03_56_B.plantData.V_body[0];
  mainV03_56_B.V_body[1] = mainV03_56_B.plantData.V_body[1];
  mainV03_56_B.V_body[2] = mainV03_56_B.plantData.V_body[2];
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
  }

  /* SignalConversion: '<S7>/SigConversion_InsertedFor_Bus Selector3_at_outport_1' */
  mainV03_56_B.V_ned[0] = mainV03_56_B.plantData.V_ned[0];
  mainV03_56_B.V_ned[1] = mainV03_56_B.plantData.V_ned[1];
  mainV03_56_B.V_ned[2] = mainV03_56_B.plantData.V_ned[2];
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
  }

  /* S-Function (saeroatmos): '<S12>/S-Function' */
  {
    /* S-Function Block: <S12>/S-Function */
    real_T *temp_table = (real_T *) &mainV03_56_DW.SFunction_temp_table[0];
    real_T *pres_table = (real_T *) &mainV03_56_DW.SFunction_pres_table[0];

    /* COESA */
    CalcAtmosCOESA( &mainV03_56_B.plantData.LLA.Altitude_m,
                   &mainV03_56_B.SFunction_o1, &mainV03_56_B.SFunction_o3,
                   &mainV03_56_B.SFunction_o4,
                   &mainV03_56_B.SFunction_o2, temp_table, pres_table, 1);
  }

  /* SignalConversion: '<S14>/TmpSignal ConversionAtSelectorInport1' */
  mainV03_56_B.TmpSignalConversionAtSelectorInport1[0] =
    mainV03_56_B.plantData.LLA.Latitude_deg;
  mainV03_56_B.TmpSignalConversionAtSelectorInport1[1] =
    mainV03_56_B.plantData.LLA.Longitude_deg;
  mainV03_56_B.TmpSignalConversionAtSelectorInport1[2] =
    mainV03_56_B.plantData.LLA.Altitude_m;

  /* Selector: '<S14>/Selector2' */
  mainV03_56_B.Selector2 = mainV03_56_B.TmpSignalConversionAtSelectorInport1[0];

  /* UnitConversion: '<S44>/Unit Conversion' */
  /* Unit Conversion - from: deg to: rad
     Expression: output = (0.0174533*input) + (0) */
  rtb_UnaryMinus_o = 0.017453292519943295 * mainV03_56_B.Selector2;

  /* Selector: '<S14>/Selector1' */
  mainV03_56_B.Selector1 = mainV03_56_B.TmpSignalConversionAtSelectorInport1[1];

  /* UnitConversion: '<S46>/Unit Conversion' */
  /* Unit Conversion - from: deg to: rad
     Expression: output = (0.0174533*input) + (0) */
  rtb_UnaryMinus_j = 0.017453292519943295 * mainV03_56_B.Selector1;

  /* Selector: '<S14>/Selector' */
  mainV03_56_B.Selector = mainV03_56_B.TmpSignalConversionAtSelectorInport1[2];

  /* Gain: '<S2>/kclock' incorporates:
   *  Clock: '<S2>/Clock'
   */
  mainV03_56_B.kclock = mainV03_56_P.kclock_Gain * mainV03_56_M->Timing.t[0];

  /* Sum: '<S2>/Sum1' incorporates:
   *  Constant: '<S2>/Initial Julian  date'
   */
  mainV03_56_B.Sum1_j = mainV03_56_P.InitialJuliandate_Value +
    mainV03_56_B.kclock;

  /* S-Function Block: <S14>/WGS84 Gravity S-Function */
  {
    int_T i;
    real_T GM, opt_m2ft, deg2rad;
    real_T *phi_ptr, *height_ptr;
    int_T uDataPrecessing, uDataCentrifugal;
    real_T uDataYear, uDataDay, uDataMonth;
    real_T E2 = WGS84_EL*WGS84_EL;     /* Linear Eccentricity Squared */
    boolean_T phi_wrapped = false;
    real_T *phi = (real_T *) &mainV03_56_DW.WGS84GravitySFunction_phi;
    real_T *h = (real_T *) &mainV03_56_DW.WGS84GravitySFunction_h;
    height_ptr = (real_T *) &mainV03_56_B.Selector;
    phi_ptr = (real_T *) &rtb_UnaryMinus_o;
    uDataYear = 99.0;
    uDataDay = mainV03_56_B.Sum1_j;

    /* get unit conversion factor  */
    opt_m2ft = 1.0;
    deg2rad = 1.0;

    /* Use Earth's Atmosphere in Gravitational Const? */
    GM = ( 1.0 == 0 ) ? WGS84_GM_DEF : WGS84_GM_PRM;

    /* WGS84EXACT */
    {
      real_T *lambda_ptr;
      real_T *gamma_h = (real_T *) &mainV03_56_DW.WGS84GravitySFunction_gamma_h;
      real_T *lambda = (real_T *) &mainV03_56_DW.WGS84GravitySFunction_lambda;
      lambda_ptr = (real_T *) &rtb_UnaryMinus_j;

      /* Get additional inputs for methods */
      for (i = 0; i < 1; i++ ) {
        real_T fphi, pi_2, flambda;

        /* create short variables for latitude (phi) and height (h) */
        phi[i] = phi_ptr[i] * deg2rad;
        h[i] = height_ptr[i] / opt_m2ft;

        /* create short variable for longitude (lambda) */
        lambda[i] = lambda_ptr[i]*deg2rad;
        if (phi[i] > AERO_PI || phi[i] < -AERO_PI) {
          phi[i] = fmod(phi[i]+AERO_PI, 2.0*AERO_PI) - AERO_PI;
        }

        fphi = fabs(phi[i]);
        pi_2 = AERO_PI/2.0;
        if (fphi > pi_2 ) {
          real_T sign = phi[i]/fphi;
          phi_wrapped = true;
          phi[i] = sign*(pi_2 - (fphi - pi_2));
        } else {
          phi_wrapped = false;
        }

        /* check and fix angle wrapping in longitude (lambda) */
        flambda = fabs(lambda[i]);
        if (phi_wrapped ) {
          lambda[i] = lambda[i] + AERO_PI;
        }

        /* check and fix angle wrapping in longitude (lambda) */
        if (flambda > AERO_PI ) {
          real_T sign = lambda[i]/flambda;
          lambda[i] = ((flambda - AERO_PI) - AERO_PI)*sign;
        }
      }

      /* Return total normal gravity */

      /* Assign udata values*/
      uDataPrecessing = (int_T) 1.0;
      uDataMonth = 10.0;
      uDataCentrifugal = (int_T) 1.0;
      wgs84_exact( h, phi, lambda, uDataMonth, uDataDay, uDataYear,
                  uDataPrecessing, uDataCentrifugal, E2, GM, opt_m2ft,
                  gamma_h, &mainV03_56_B.MatrixConcatenate[2],
                  &mainV03_56_B.MatrixConcatenate[0], 1);

      /* east component is zero */
      mainV03_56_B.MatrixConcatenate[1] =
        mainV03_56_DW.WGS84GravitySFunction_gamma_phi;
    }
  }

  /* Sum: '<S2>/Sum' incorporates:
   *  Constant: '<S2>/Initial Decimal date'
   */
  mainV03_56_B.Sum_p = mainV03_56_B.kclock +
    mainV03_56_P.InitialDecimaldate_Value;
  for (s155_iter = 0; s155_iter < 3; s155_iter++) {
    /* Gain: '<S2>/Gravity in Earth Axes' */
    mainV03_56_B.GravityinEarthAxes[s155_iter] =
      mainV03_56_P.GravityinEarthAxes_Gain[s155_iter] *
      mainV03_56_B.MatrixConcatenate[s155_iter];

    /* Math: '<S2>/Math Function' */
    mainV03_56_B.MathFunction[s155_iter] =
      mainV03_56_B.GravityinEarthAxes[s155_iter];

    /* Sum: '<S13>/Sum2' incorporates:
     *  Constant: '<S13>/motor2 position '
     */
    mainV03_56_B.r_A_body_1[s155_iter] =
      mainV03_56_P.motor2position_Value_c[s155_iter] -
      mainV03_56_B.plantData.CG[s155_iter];

    /* Math: '<S13>/Transpose' */
    mainV03_56_B.Transpose_m[3 * s155_iter] =
      mainV03_56_B.plantData.DCM_body_earth[s155_iter];
    mainV03_56_B.Transpose_m[1 + 3 * s155_iter] =
      mainV03_56_B.plantData.DCM_body_earth[s155_iter + 3];
    mainV03_56_B.Transpose_m[2 + 3 * s155_iter] =
      mainV03_56_B.plantData.DCM_body_earth[s155_iter + 6];
  }

  /* Product: '<S31>/i x j' */
  mainV03_56_B.ixj = mainV03_56_B.plantData.Omega_body.p *
    mainV03_56_B.r_A_body_1[1];

  /* Product: '<S31>/j x k' */
  mainV03_56_B.jxk = mainV03_56_B.plantData.Omega_body.q *
    mainV03_56_B.r_A_body_1[2];

  /* Product: '<S31>/k x i' */
  mainV03_56_B.kxi = mainV03_56_B.plantData.Omega_body.r *
    mainV03_56_B.r_A_body_1[0];

  /* Product: '<S32>/i x k' */
  mainV03_56_B.ixk = mainV03_56_B.plantData.Omega_body.p *
    mainV03_56_B.r_A_body_1[2];

  /* Product: '<S32>/j x i' */
  mainV03_56_B.jxi = mainV03_56_B.plantData.Omega_body.q *
    mainV03_56_B.r_A_body_1[0];

  /* Product: '<S32>/k x j' */
  mainV03_56_B.kxj = mainV03_56_B.plantData.Omega_body.r *
    mainV03_56_B.r_A_body_1[1];

  /* Sum: '<S22>/Sum' */
  mainV03_56_B.Sum_e[0] = mainV03_56_B.jxk - mainV03_56_B.kxj;
  mainV03_56_B.Sum_e[1] = mainV03_56_B.kxi - mainV03_56_B.ixk;
  mainV03_56_B.Sum_e[2] = mainV03_56_B.ixj - mainV03_56_B.jxi;
  for (s155_iter = 0; s155_iter < 3; s155_iter++) {
    /* Product: '<S13>/To Earth Axes2' */
    mainV03_56_B.r_A_earth_1[s155_iter] = 0.0;
    mainV03_56_B.r_A_earth_1[s155_iter] += mainV03_56_B.Transpose_m[s155_iter] *
      mainV03_56_B.r_A_body_1[0];
    mainV03_56_B.r_A_earth_1[s155_iter] += mainV03_56_B.Transpose_m[s155_iter +
      3] * mainV03_56_B.r_A_body_1[1];
    mainV03_56_B.r_A_earth_1[s155_iter] += mainV03_56_B.Transpose_m[s155_iter +
      6] * mainV03_56_B.r_A_body_1[2];

    /* Sum: '<S13>/Sum3' */
    mainV03_56_B.r_A_earth_0[s155_iter] = mainV03_56_B.r_A_earth_1[s155_iter] +
      mainV03_56_B.plantData.X_ned[s155_iter];

    /* Product: '<S13>/To Earth Axes1' */
    mainV03_56_B.v_A_CG_0[s155_iter] = 0.0;
    mainV03_56_B.v_A_CG_0[s155_iter] += mainV03_56_B.Transpose_m[s155_iter] *
      mainV03_56_B.Sum_e[0];
    mainV03_56_B.v_A_CG_0[s155_iter] += mainV03_56_B.Transpose_m[s155_iter + 3] *
      mainV03_56_B.Sum_e[1];
    mainV03_56_B.v_A_CG_0[s155_iter] += mainV03_56_B.Transpose_m[s155_iter + 6] *
      mainV03_56_B.Sum_e[2];

    /* Sum: '<S13>/Sum4' */
    mainV03_56_B.v_A_earth_0[s155_iter] = mainV03_56_B.v_A_CG_0[s155_iter] +
      mainV03_56_B.plantData.V_ned[s155_iter];
  }

  /* Sum: '<S28>/Sum5' incorporates:
   *  Constant: '<S28>/hground_value1'
   *  UnaryMinus: '<S28>/Ze2height1'
   */
  mainV03_56_B.z = -mainV03_56_B.r_A_earth_0[2] -
    mainV03_56_P.hground_value1_Value;

  /* Switch: '<S28>/Switch3' incorporates:
   *  Constant: '<S28>/Constant3'
   */
  if (mainV03_56_B.z > mainV03_56_P.Switch3_Threshold) {
    mainV03_56_B.Switch3 = mainV03_56_P.Constant3_Value_lv;
  } else {
    mainV03_56_B.Switch3 = mainV03_56_B.z;
  }

  /* End of Switch: '<S28>/Switch3' */

  /* UnaryMinus: '<S28>/height2d2' */
  rtb_UnaryMinus_j = -mainV03_56_B.Switch3;

  /* Sqrt: '<S28>/Sqrt1' */
  rtb_UnaryMinus_j = sqrt(rtb_UnaryMinus_j);

  /* Gain: '<S28>/-K2' */
  mainV03_56_B.K2 = -mainV03_56_P.groundReaction.K * mainV03_56_B.z;

  /* Gain: '<S28>/-D1' incorporates:
   *  UnaryMinus: '<S28>/Ze_dot2height_dot3'
   */
  mainV03_56_B.D1 = -mainV03_56_P.groundReaction.D * -mainV03_56_B.v_A_earth_0[2];

  /* Sum: '<S28>/Sum4' */
  mainV03_56_B.Sum4 = mainV03_56_B.K2 + mainV03_56_B.D1;

  /* Product: '<S28>/Product1' */
  mainV03_56_B.Product1_a = rtb_UnaryMinus_j * mainV03_56_B.Sum4;

  /* UnaryMinus: '<S28>/height2d3' */
  rtb_UnaryMinus_j = -mainV03_56_B.Product1_a;

  /* SignalConversion: '<S13>/TmpSignal ConversionAtTo Body Axes3Inport2' incorporates:
   *  Constant: '<S28>/Constant4'
   *  Constant: '<S28>/Constant5'
   */
  mainV03_56_B.TmpSignalConversionAtToBodyAxes3Inport2[0] =
    mainV03_56_P.Constant4_Value_p;
  mainV03_56_B.TmpSignalConversionAtToBodyAxes3Inport2[1] =
    mainV03_56_P.Constant5_Value_l;
  mainV03_56_B.TmpSignalConversionAtToBodyAxes3Inport2[2] = rtb_UnaryMinus_j;
  for (s155_iter = 0; s155_iter < 3; s155_iter++) {
    /* Product: '<S13>/To Body Axes3' */
    mainV03_56_B.F_A_body[s155_iter] = 0.0;
    mainV03_56_B.F_A_body[s155_iter] +=
      mainV03_56_B.plantData.DCM_body_earth[s155_iter] *
      mainV03_56_B.TmpSignalConversionAtToBodyAxes3Inport2[0];
    mainV03_56_B.F_A_body[s155_iter] +=
      mainV03_56_B.plantData.DCM_body_earth[s155_iter + 3] *
      mainV03_56_B.TmpSignalConversionAtToBodyAxes3Inport2[1];
    mainV03_56_B.F_A_body[s155_iter] +=
      mainV03_56_B.plantData.DCM_body_earth[s155_iter + 6] *
      mainV03_56_B.TmpSignalConversionAtToBodyAxes3Inport2[2];

    /* Sum: '<S13>/Sum6' incorporates:
     *  Constant: '<S13>/motor2 position 1'
     */
    mainV03_56_B.r_A_body_1_g[s155_iter] =
      mainV03_56_P.motor2position1_Value_m[s155_iter] -
      mainV03_56_B.plantData.CG[s155_iter];
  }

  /* Product: '<S33>/i x j' */
  mainV03_56_B.ixj_j = mainV03_56_B.r_A_body_1[0] * mainV03_56_B.F_A_body[1];

  /* Product: '<S33>/j x k' */
  mainV03_56_B.jxk_d = mainV03_56_B.r_A_body_1[1] * mainV03_56_B.F_A_body[2];

  /* Product: '<S33>/k x i' */
  mainV03_56_B.kxi_c = mainV03_56_B.r_A_body_1[2] * mainV03_56_B.F_A_body[0];

  /* Product: '<S34>/i x k' */
  mainV03_56_B.ixk_n = mainV03_56_B.r_A_body_1[0] * mainV03_56_B.F_A_body[2];

  /* Product: '<S34>/j x i' */
  mainV03_56_B.jxi_g = mainV03_56_B.r_A_body_1[1] * mainV03_56_B.F_A_body[0];

  /* Product: '<S34>/k x j' */
  mainV03_56_B.kxj_m = mainV03_56_B.r_A_body_1[2] * mainV03_56_B.F_A_body[1];

  /* Sum: '<S23>/Sum' */
  mainV03_56_B.Sum_d[0] = mainV03_56_B.jxk_d - mainV03_56_B.kxj_m;
  mainV03_56_B.Sum_d[1] = mainV03_56_B.kxi_c - mainV03_56_B.ixk_n;
  mainV03_56_B.Sum_d[2] = mainV03_56_B.ixj_j - mainV03_56_B.jxi_g;

  /* Product: '<S35>/i x j' */
  mainV03_56_B.ixj_i = mainV03_56_B.plantData.Omega_body.p *
    mainV03_56_B.r_A_body_1_g[1];

  /* Product: '<S35>/j x k' */
  mainV03_56_B.jxk_p = mainV03_56_B.plantData.Omega_body.q *
    mainV03_56_B.r_A_body_1_g[2];

  /* Product: '<S35>/k x i' */
  mainV03_56_B.kxi_i = mainV03_56_B.plantData.Omega_body.r *
    mainV03_56_B.r_A_body_1_g[0];

  /* Product: '<S36>/i x k' */
  mainV03_56_B.ixk_o = mainV03_56_B.plantData.Omega_body.p *
    mainV03_56_B.r_A_body_1_g[2];

  /* Product: '<S36>/j x i' */
  mainV03_56_B.jxi_f = mainV03_56_B.plantData.Omega_body.q *
    mainV03_56_B.r_A_body_1_g[0];

  /* Product: '<S36>/k x j' */
  mainV03_56_B.kxj_i = mainV03_56_B.plantData.Omega_body.r *
    mainV03_56_B.r_A_body_1_g[1];

  /* Sum: '<S24>/Sum' */
  mainV03_56_B.Sum_k[0] = mainV03_56_B.jxk_p - mainV03_56_B.kxj_i;
  mainV03_56_B.Sum_k[1] = mainV03_56_B.kxi_i - mainV03_56_B.ixk_o;
  mainV03_56_B.Sum_k[2] = mainV03_56_B.ixj_i - mainV03_56_B.jxi_f;
  for (s155_iter = 0; s155_iter < 3; s155_iter++) {
    /* Product: '<S13>/To Earth Axes4' */
    mainV03_56_B.r_A_earth_1_b[s155_iter] = 0.0;
    mainV03_56_B.r_A_earth_1_b[s155_iter] += mainV03_56_B.Transpose_m[s155_iter]
      * mainV03_56_B.r_A_body_1_g[0];
    mainV03_56_B.r_A_earth_1_b[s155_iter] += mainV03_56_B.Transpose_m[s155_iter
      + 3] * mainV03_56_B.r_A_body_1_g[1];
    mainV03_56_B.r_A_earth_1_b[s155_iter] += mainV03_56_B.Transpose_m[s155_iter
      + 6] * mainV03_56_B.r_A_body_1_g[2];

    /* Sum: '<S13>/Sum7' */
    mainV03_56_B.r_A_earth_0_n[s155_iter] = mainV03_56_B.r_A_earth_1_b[s155_iter]
      + mainV03_56_B.plantData.X_ned[s155_iter];

    /* Product: '<S13>/To Earth Axes3' */
    mainV03_56_B.v_A_CG_0_l[s155_iter] = 0.0;
    mainV03_56_B.v_A_CG_0_l[s155_iter] += mainV03_56_B.Transpose_m[s155_iter] *
      mainV03_56_B.Sum_k[0];
    mainV03_56_B.v_A_CG_0_l[s155_iter] += mainV03_56_B.Transpose_m[s155_iter + 3]
      * mainV03_56_B.Sum_k[1];
    mainV03_56_B.v_A_CG_0_l[s155_iter] += mainV03_56_B.Transpose_m[s155_iter + 6]
      * mainV03_56_B.Sum_k[2];

    /* Sum: '<S13>/Sum8' */
    mainV03_56_B.v_A_earth_0_b[s155_iter] = mainV03_56_B.v_A_CG_0_l[s155_iter] +
      mainV03_56_B.plantData.V_ned[s155_iter];
  }

  /* Sum: '<S29>/Sum5' incorporates:
   *  Constant: '<S29>/hground_value1'
   *  UnaryMinus: '<S29>/Ze2height1'
   */
  mainV03_56_B.z_l = -mainV03_56_B.r_A_earth_0_n[2] -
    mainV03_56_P.hground_value1_Value_o;

  /* Switch: '<S29>/Switch3' incorporates:
   *  Constant: '<S29>/Constant3'
   */
  if (mainV03_56_B.z_l > mainV03_56_P.Switch3_Threshold_g) {
    mainV03_56_B.Switch3_k = mainV03_56_P.Constant3_Value_k;
  } else {
    mainV03_56_B.Switch3_k = mainV03_56_B.z_l;
  }

  /* End of Switch: '<S29>/Switch3' */

  /* UnaryMinus: '<S29>/height2d2' */
  rtb_UnaryMinus_j = -mainV03_56_B.Switch3_k;

  /* Sqrt: '<S29>/Sqrt1' */
  rtb_UnaryMinus_j = sqrt(rtb_UnaryMinus_j);

  /* Gain: '<S29>/-K2' */
  mainV03_56_B.K2_n = -mainV03_56_P.groundReaction.K * mainV03_56_B.z_l;

  /* Gain: '<S29>/-D1' incorporates:
   *  UnaryMinus: '<S29>/Ze_dot2height_dot3'
   */
  mainV03_56_B.D1_e = -mainV03_56_P.groundReaction.D *
    -mainV03_56_B.v_A_earth_0_b[2];

  /* Sum: '<S29>/Sum4' */
  mainV03_56_B.Sum4_j = mainV03_56_B.K2_n + mainV03_56_B.D1_e;

  /* Product: '<S29>/Product1' */
  mainV03_56_B.Product1_d = rtb_UnaryMinus_j * mainV03_56_B.Sum4_j;

  /* UnaryMinus: '<S29>/height2d3' */
  rtb_UnaryMinus_j = -mainV03_56_B.Product1_d;

  /* SignalConversion: '<S13>/TmpSignal ConversionAtTo Body Axes1Inport2' incorporates:
   *  Constant: '<S29>/Constant4'
   *  Constant: '<S29>/Constant5'
   */
  mainV03_56_B.TmpSignalConversionAtToBodyAxes1Inport2[0] =
    mainV03_56_P.Constant4_Value_c;
  mainV03_56_B.TmpSignalConversionAtToBodyAxes1Inport2[1] =
    mainV03_56_P.Constant5_Value_b;
  mainV03_56_B.TmpSignalConversionAtToBodyAxes1Inport2[2] = rtb_UnaryMinus_j;
  for (s155_iter = 0; s155_iter < 3; s155_iter++) {
    /* Product: '<S13>/To Body Axes1' */
    mainV03_56_B.F_B_body[s155_iter] = 0.0;
    mainV03_56_B.F_B_body[s155_iter] +=
      mainV03_56_B.plantData.DCM_body_earth[s155_iter] *
      mainV03_56_B.TmpSignalConversionAtToBodyAxes1Inport2[0];
    mainV03_56_B.F_B_body[s155_iter] +=
      mainV03_56_B.plantData.DCM_body_earth[s155_iter + 3] *
      mainV03_56_B.TmpSignalConversionAtToBodyAxes1Inport2[1];
    mainV03_56_B.F_B_body[s155_iter] +=
      mainV03_56_B.plantData.DCM_body_earth[s155_iter + 6] *
      mainV03_56_B.TmpSignalConversionAtToBodyAxes1Inport2[2];

    /* Sum: '<S13>/Sum9' incorporates:
     *  Constant: '<S13>/motor2 position 2'
     */
    mainV03_56_B.r_A_body_1_e[s155_iter] =
      mainV03_56_P.motor2position2_Value[s155_iter] -
      mainV03_56_B.plantData.CG[s155_iter];
  }

  /* Product: '<S37>/i x j' */
  mainV03_56_B.ixj_m = mainV03_56_B.r_A_body_1_g[0] * mainV03_56_B.F_B_body[1];

  /* Product: '<S37>/j x k' */
  mainV03_56_B.jxk_k = mainV03_56_B.r_A_body_1_g[1] * mainV03_56_B.F_B_body[2];

  /* Product: '<S37>/k x i' */
  mainV03_56_B.kxi_l = mainV03_56_B.r_A_body_1_g[2] * mainV03_56_B.F_B_body[0];

  /* Product: '<S38>/i x k' */
  mainV03_56_B.ixk_d = mainV03_56_B.r_A_body_1_g[0] * mainV03_56_B.F_B_body[2];

  /* Product: '<S38>/j x i' */
  mainV03_56_B.jxi_k = mainV03_56_B.r_A_body_1_g[1] * mainV03_56_B.F_B_body[0];

  /* Product: '<S38>/k x j' */
  mainV03_56_B.kxj_mr = mainV03_56_B.r_A_body_1_g[2] * mainV03_56_B.F_B_body[1];

  /* Sum: '<S25>/Sum' */
  mainV03_56_B.Sum_a[0] = mainV03_56_B.jxk_k - mainV03_56_B.kxj_mr;
  mainV03_56_B.Sum_a[1] = mainV03_56_B.kxi_l - mainV03_56_B.ixk_d;
  mainV03_56_B.Sum_a[2] = mainV03_56_B.ixj_m - mainV03_56_B.jxi_k;

  /* Product: '<S39>/i x j' */
  mainV03_56_B.ixj_b = mainV03_56_B.plantData.Omega_body.p *
    mainV03_56_B.r_A_body_1_e[1];

  /* Product: '<S39>/j x k' */
  mainV03_56_B.jxk_c = mainV03_56_B.plantData.Omega_body.q *
    mainV03_56_B.r_A_body_1_e[2];

  /* Product: '<S39>/k x i' */
  mainV03_56_B.kxi_g = mainV03_56_B.plantData.Omega_body.r *
    mainV03_56_B.r_A_body_1_e[0];

  /* Product: '<S40>/i x k' */
  mainV03_56_B.ixk_c = mainV03_56_B.plantData.Omega_body.p *
    mainV03_56_B.r_A_body_1_e[2];

  /* Product: '<S40>/j x i' */
  mainV03_56_B.jxi_j = mainV03_56_B.plantData.Omega_body.q *
    mainV03_56_B.r_A_body_1_e[0];

  /* Product: '<S40>/k x j' */
  mainV03_56_B.kxj_a = mainV03_56_B.plantData.Omega_body.r *
    mainV03_56_B.r_A_body_1_e[1];

  /* Sum: '<S26>/Sum' */
  mainV03_56_B.Sum_cs[0] = mainV03_56_B.jxk_c - mainV03_56_B.kxj_a;
  mainV03_56_B.Sum_cs[1] = mainV03_56_B.kxi_g - mainV03_56_B.ixk_c;
  mainV03_56_B.Sum_cs[2] = mainV03_56_B.ixj_b - mainV03_56_B.jxi_j;
  for (s155_iter = 0; s155_iter < 3; s155_iter++) {
    /* Product: '<S13>/To Earth Axes6' */
    mainV03_56_B.r_A_earth_1_i[s155_iter] = 0.0;
    mainV03_56_B.r_A_earth_1_i[s155_iter] += mainV03_56_B.Transpose_m[s155_iter]
      * mainV03_56_B.r_A_body_1_e[0];
    mainV03_56_B.r_A_earth_1_i[s155_iter] += mainV03_56_B.Transpose_m[s155_iter
      + 3] * mainV03_56_B.r_A_body_1_e[1];
    mainV03_56_B.r_A_earth_1_i[s155_iter] += mainV03_56_B.Transpose_m[s155_iter
      + 6] * mainV03_56_B.r_A_body_1_e[2];

    /* Sum: '<S13>/Sum10' */
    mainV03_56_B.r_A_earth_0_e[s155_iter] = mainV03_56_B.r_A_earth_1_i[s155_iter]
      + mainV03_56_B.plantData.X_ned[s155_iter];

    /* Product: '<S13>/To Earth Axes5' */
    mainV03_56_B.v_A_CG_0_g[s155_iter] = 0.0;
    mainV03_56_B.v_A_CG_0_g[s155_iter] += mainV03_56_B.Transpose_m[s155_iter] *
      mainV03_56_B.Sum_cs[0];
    mainV03_56_B.v_A_CG_0_g[s155_iter] += mainV03_56_B.Transpose_m[s155_iter + 3]
      * mainV03_56_B.Sum_cs[1];
    mainV03_56_B.v_A_CG_0_g[s155_iter] += mainV03_56_B.Transpose_m[s155_iter + 6]
      * mainV03_56_B.Sum_cs[2];

    /* Sum: '<S13>/Sum11' */
    mainV03_56_B.v_A_earth_0_h[s155_iter] = mainV03_56_B.v_A_CG_0_g[s155_iter] +
      mainV03_56_B.plantData.V_ned[s155_iter];
  }

  /* Sum: '<S30>/Sum5' incorporates:
   *  Constant: '<S30>/hground_value1'
   *  UnaryMinus: '<S30>/Ze2height1'
   */
  mainV03_56_B.z_n = -mainV03_56_B.r_A_earth_0_e[2] -
    mainV03_56_P.hground_value1_Value_h;

  /* Switch: '<S30>/Switch3' incorporates:
   *  Constant: '<S30>/Constant3'
   */
  if (mainV03_56_B.z_n > mainV03_56_P.Switch3_Threshold_b) {
    mainV03_56_B.Switch3_ks = mainV03_56_P.Constant3_Value_h;
  } else {
    mainV03_56_B.Switch3_ks = mainV03_56_B.z_n;
  }

  /* End of Switch: '<S30>/Switch3' */

  /* UnaryMinus: '<S30>/height2d2' */
  rtb_UnaryMinus_j = -mainV03_56_B.Switch3_ks;

  /* Sqrt: '<S30>/Sqrt1' */
  rtb_UnaryMinus_j = sqrt(rtb_UnaryMinus_j);

  /* Gain: '<S30>/-K2' */
  mainV03_56_B.K2_l = -mainV03_56_P.groundReaction.K * mainV03_56_B.z_n;

  /* Gain: '<S30>/-D1' incorporates:
   *  UnaryMinus: '<S30>/Ze_dot2height_dot3'
   */
  mainV03_56_B.D1_o = -mainV03_56_P.groundReaction.D *
    -mainV03_56_B.v_A_earth_0_h[2];

  /* Sum: '<S30>/Sum4' */
  mainV03_56_B.Sum4_f = mainV03_56_B.K2_l + mainV03_56_B.D1_o;

  /* Product: '<S30>/Product1' */
  mainV03_56_B.Product1_g = rtb_UnaryMinus_j * mainV03_56_B.Sum4_f;

  /* UnaryMinus: '<S30>/height2d3' */
  rtb_UnaryMinus_j = -mainV03_56_B.Product1_g;

  /* SignalConversion: '<S13>/TmpSignal ConversionAtTo Body Axes2Inport2' incorporates:
   *  Constant: '<S30>/Constant4'
   *  Constant: '<S30>/Constant5'
   */
  mainV03_56_B.TmpSignalConversionAtToBodyAxes2Inport2[0] =
    mainV03_56_P.Constant4_Value_c3;
  mainV03_56_B.TmpSignalConversionAtToBodyAxes2Inport2[1] =
    mainV03_56_P.Constant5_Value_d;
  mainV03_56_B.TmpSignalConversionAtToBodyAxes2Inport2[2] = rtb_UnaryMinus_j;

  /* Product: '<S13>/To Body Axes2' */
  for (rtb_Sum1_eh = 0; rtb_Sum1_eh < 3; rtb_Sum1_eh++) {
    mainV03_56_B.F_C_body[rtb_Sum1_eh] = 0.0;
    mainV03_56_B.F_C_body[rtb_Sum1_eh] +=
      mainV03_56_B.plantData.DCM_body_earth[rtb_Sum1_eh] *
      mainV03_56_B.TmpSignalConversionAtToBodyAxes2Inport2[0];
    mainV03_56_B.F_C_body[rtb_Sum1_eh] +=
      mainV03_56_B.plantData.DCM_body_earth[rtb_Sum1_eh + 3] *
      mainV03_56_B.TmpSignalConversionAtToBodyAxes2Inport2[1];
    mainV03_56_B.F_C_body[rtb_Sum1_eh] +=
      mainV03_56_B.plantData.DCM_body_earth[rtb_Sum1_eh + 6] *
      mainV03_56_B.TmpSignalConversionAtToBodyAxes2Inport2[2];
  }

  /* End of Product: '<S13>/To Body Axes2' */

  /* Product: '<S41>/i x j' */
  mainV03_56_B.ixj_g = mainV03_56_B.r_A_body_1_e[0] * mainV03_56_B.F_C_body[1];

  /* Product: '<S41>/j x k' */
  mainV03_56_B.jxk_i = mainV03_56_B.r_A_body_1_e[1] * mainV03_56_B.F_C_body[2];

  /* Product: '<S41>/k x i' */
  mainV03_56_B.kxi_a = mainV03_56_B.r_A_body_1_e[2] * mainV03_56_B.F_C_body[0];

  /* Product: '<S42>/i x k' */
  mainV03_56_B.ixk_f = mainV03_56_B.r_A_body_1_e[0] * mainV03_56_B.F_C_body[2];

  /* Product: '<S42>/j x i' */
  mainV03_56_B.jxi_b = mainV03_56_B.r_A_body_1_e[1] * mainV03_56_B.F_C_body[0];

  /* Product: '<S42>/k x j' */
  mainV03_56_B.kxj_b = mainV03_56_B.r_A_body_1_e[2] * mainV03_56_B.F_C_body[1];

  /* Sum: '<S27>/Sum' */
  mainV03_56_B.Sum_f[0] = mainV03_56_B.jxk_i - mainV03_56_B.kxj_b;
  mainV03_56_B.Sum_f[1] = mainV03_56_B.kxi_a - mainV03_56_B.ixk_f;
  mainV03_56_B.Sum_f[2] = mainV03_56_B.ixj_g - mainV03_56_B.jxi_b;

  /* InitialCondition: '<S13>/IC' */
  if ((mainV03_56_DW.IC_FirstOutputTime == (rtMinusInf)) ||
      (mainV03_56_DW.IC_FirstOutputTime == mainV03_56_M->Timing.t[0])) {
    mainV03_56_DW.IC_FirstOutputTime = mainV03_56_M->Timing.t[0];
    rtb_sincos_o2_p[0] = mainV03_56_P.IC_Value_e[0];
    rtb_sincos_o2_p[1] = mainV03_56_P.IC_Value_e[1];
    rtb_sincos_o2_p[2] = mainV03_56_P.IC_Value_e[2];
  } else {
    rtb_sincos_o2_p[0] = mainV03_56_B.Sum_d[0];
    rtb_sincos_o2_p[1] = mainV03_56_B.Sum_d[1];
    rtb_sincos_o2_p[2] = mainV03_56_B.Sum_d[2];
  }

  /* End of InitialCondition: '<S13>/IC' */

  /* InitialCondition: '<S13>/IC1' */
  if ((mainV03_56_DW.IC1_FirstOutputTime == (rtMinusInf)) ||
      (mainV03_56_DW.IC1_FirstOutputTime == mainV03_56_M->Timing.t[0])) {
    mainV03_56_DW.IC1_FirstOutputTime = mainV03_56_M->Timing.t[0];
    rtb_sincos_o1_n_idx_0 = mainV03_56_P.IC1_Value_b[0];
    rtb_sincos_o1_n_idx_1 = mainV03_56_P.IC1_Value_b[1];
    rtb_Sum1_no = mainV03_56_P.IC1_Value_b[2];
  } else {
    rtb_sincos_o1_n_idx_0 = mainV03_56_B.F_A_body[0];
    rtb_sincos_o1_n_idx_1 = mainV03_56_B.F_A_body[1];
    rtb_Sum1_no = mainV03_56_B.F_A_body[2];
  }

  /* End of InitialCondition: '<S13>/IC1' */

  /* InitialCondition: '<S13>/IC2' */
  if ((mainV03_56_DW.IC2_FirstOutputTime_d == (rtMinusInf)) ||
      (mainV03_56_DW.IC2_FirstOutputTime_d == mainV03_56_M->Timing.t[0])) {
    mainV03_56_DW.IC2_FirstOutputTime_d = mainV03_56_M->Timing.t[0];
    rtb_sincos_o2_i_idx_0 = mainV03_56_P.IC2_Value_c[0];
    rtb_sincos_o2_i_idx_1 = mainV03_56_P.IC2_Value_c[1];
    rtb_sincos_o2_i_idx_2 = mainV03_56_P.IC2_Value_c[2];
  } else {
    rtb_sincos_o2_i_idx_0 = mainV03_56_B.Sum_a[0];
    rtb_sincos_o2_i_idx_1 = mainV03_56_B.Sum_a[1];
    rtb_sincos_o2_i_idx_2 = mainV03_56_B.Sum_a[2];
  }

  /* End of InitialCondition: '<S13>/IC2' */

  /* InitialCondition: '<S13>/IC3' */
  if ((mainV03_56_DW.IC3_FirstOutputTime_o == (rtMinusInf)) ||
      (mainV03_56_DW.IC3_FirstOutputTime_o == mainV03_56_M->Timing.t[0])) {
    mainV03_56_DW.IC3_FirstOutputTime_o = mainV03_56_M->Timing.t[0];
    rtb_IC3_idx_0 = mainV03_56_P.IC3_Value_k[0];
    rtb_IC3_idx_1 = mainV03_56_P.IC3_Value_k[1];
    rtb_IC3_idx_2 = mainV03_56_P.IC3_Value_k[2];
  } else {
    rtb_IC3_idx_0 = mainV03_56_B.F_B_body[0];
    rtb_IC3_idx_1 = mainV03_56_B.F_B_body[1];
    rtb_IC3_idx_2 = mainV03_56_B.F_B_body[2];
  }

  /* End of InitialCondition: '<S13>/IC3' */

  /* InitialCondition: '<S13>/IC4' */
  if ((mainV03_56_DW.IC4_FirstOutputTime_m == (rtMinusInf)) ||
      (mainV03_56_DW.IC4_FirstOutputTime_m == mainV03_56_M->Timing.t[0])) {
    mainV03_56_DW.IC4_FirstOutputTime_m = mainV03_56_M->Timing.t[0];
    rtb_IC4_idx_0 = mainV03_56_P.IC4_Value_g[0];
    rtb_IC4_idx_1 = mainV03_56_P.IC4_Value_g[1];
    rtb_IC4_idx_2 = mainV03_56_P.IC4_Value_g[2];
  } else {
    rtb_IC4_idx_0 = mainV03_56_B.Sum_f[0];
    rtb_IC4_idx_1 = mainV03_56_B.Sum_f[1];
    rtb_IC4_idx_2 = mainV03_56_B.Sum_f[2];
  }

  /* End of InitialCondition: '<S13>/IC4' */

  /* InitialCondition: '<S13>/IC5' */
  if ((mainV03_56_DW.IC5_FirstOutputTime_d == (rtMinusInf)) ||
      (mainV03_56_DW.IC5_FirstOutputTime_d == mainV03_56_M->Timing.t[0])) {
    mainV03_56_DW.IC5_FirstOutputTime_d = mainV03_56_M->Timing.t[0];
    rtb_Shiftright[0] = mainV03_56_P.IC5_Value_b[0];
    rtb_Shiftright[1] = mainV03_56_P.IC5_Value_b[1];
    rtb_Shiftright[2] = mainV03_56_P.IC5_Value_b[2];
  } else {
    rtb_Shiftright[0] = mainV03_56_B.F_C_body[0];
    rtb_Shiftright[1] = mainV03_56_B.F_C_body[1];
    rtb_Shiftright[2] = mainV03_56_B.F_C_body[2];
  }

  /* End of InitialCondition: '<S13>/IC5' */

  /* Sum: '<S13>/Sum1' */
  mainV03_56_B.envData.Fground[0] = (rtb_sincos_o1_n_idx_0 + rtb_IC3_idx_0) +
    rtb_Shiftright[0];

  /* Sum: '<S13>/Sum5' */
  mainV03_56_B.envData.Mground[0] = (rtb_sincos_o2_p[0] + rtb_sincos_o2_i_idx_0)
    + rtb_IC4_idx_0;

  /* Sum: '<S13>/Sum1' */
  mainV03_56_B.envData.Fground[1] = (rtb_sincos_o1_n_idx_1 + rtb_IC3_idx_1) +
    rtb_Shiftright[1];

  /* Sum: '<S13>/Sum5' */
  mainV03_56_B.envData.Mground[1] = (rtb_sincos_o2_p[1] + rtb_sincos_o2_i_idx_1)
    + rtb_IC4_idx_1;

  /* Sum: '<S13>/Sum1' */
  mainV03_56_B.envData.Fground[2] = (rtb_Sum1_no + rtb_IC3_idx_2) +
    rtb_Shiftright[2];

  /* Sum: '<S13>/Sum5' */
  mainV03_56_B.envData.Mground[2] = (rtb_sincos_o2_p[2] + rtb_sincos_o2_i_idx_2)
    + rtb_IC4_idx_2;

  /* Sum: '<S104>/Sum' incorporates:
   *  Constant: '<S104>/epoch'
   */
  mainV03_56_B.Sum_n = mainV03_56_B.Sum_p - mainV03_56_P.epoch_Value;
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Memory: '<S112>/otime' */
    mainV03_56_B.otime = mainV03_56_DW.otime_PreviousInput;
  }

  /* RelationalOperator: '<S112>/Relational Operator' */
  mainV03_56_B.RelationalOperator = (mainV03_56_B.Sum_p != mainV03_56_B.otime);

  /* Saturate: '<S16>/+//- 180 deg' */
  if (mainV03_56_B.plantData.LLA.Longitude_deg > mainV03_56_P.u80deg_UpperSat) {
    mainV03_56_B.u80deg = mainV03_56_P.u80deg_UpperSat;
  } else if (mainV03_56_B.plantData.LLA.Longitude_deg <
             mainV03_56_P.u80deg_LowerSat) {
    mainV03_56_B.u80deg = mainV03_56_P.u80deg_LowerSat;
  } else {
    mainV03_56_B.u80deg = mainV03_56_B.plantData.LLA.Longitude_deg;
  }

  /* End of Saturate: '<S16>/+//- 180 deg' */

  /* Saturate: '<S16>/+//- 90 deg' */
  if (mainV03_56_B.plantData.LLA.Latitude_deg > mainV03_56_P.u0deg_UpperSat) {
    mainV03_56_B.u0deg = mainV03_56_P.u0deg_UpperSat;
  } else if (mainV03_56_B.plantData.LLA.Latitude_deg <
             mainV03_56_P.u0deg_LowerSat) {
    mainV03_56_B.u0deg = mainV03_56_P.u0deg_LowerSat;
  } else {
    mainV03_56_B.u0deg = mainV03_56_B.plantData.LLA.Latitude_deg;
  }

  /* End of Saturate: '<S16>/+//- 90 deg' */

  /* UnitConversion: '<S156>/Unit Conversion' */
  /* Unit Conversion - from: deg to: rad
     Expression: output = (0.0174533*input) + (0) */
  rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[0] = 0.017453292519943295
    * mainV03_56_B.u80deg;
  rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[1] = 0.017453292519943295
    * mainV03_56_B.u0deg;

  /* Trigonometry: '<S109>/sincos' */
  mainV03_56_B.sincos_o1[0] = sin
    (rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[0]);
  mainV03_56_B.sincos_o2[0] = cos
    (rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[0]);
  mainV03_56_B.sincos_o1[1] = sin
    (rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[1]);
  mainV03_56_B.sincos_o2[1] = cos
    (rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[1]);
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Memory: '<S111>/olon' */
    mainV03_56_B.olon = mainV03_56_DW.olon_PreviousInput;
  }

  /* RelationalOperator: '<S111>/Relational Operator' */
  mainV03_56_B.RelationalOperator_b = (mainV03_56_B.u80deg != mainV03_56_B.olon);
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* SignalConversion: '<S104>/HiddenBuf_InsertedFor_Convert from geodetic to  spherical coordinates _at_inport_2' */
    mainV03_56_B.HiddenBuf_InsertedFor_Convertfromgeodetictosphericalcoordinates_at_inport_2
      = mainV03_56_B.RelationalOperator_b;

    /* Outputs for Enabled SubSystem: '<S104>/Convert from geodetic to  spherical coordinates ' incorporates:
     *  EnablePort: '<S108>/Enable'
     */
    if (rtmIsMajorTimeStep(mainV03_56_M)) {
      if (mainV03_56_B.HiddenBuf_InsertedFor_Convertfromgeodetictosphericalcoordinates_at_inport_2)
      {
        if (!mainV03_56_DW.Convertfromgeodetictosphericalcoordinates_MODE) {
          mainV03_56_DW.Convertfromgeodetictosphericalcoordinates_MODE = true;
        }
      } else {
        if (mainV03_56_DW.Convertfromgeodetictosphericalcoordinates_MODE) {
          mainV03_56_DW.Convertfromgeodetictosphericalcoordinates_MODE = false;
        }
      }
    }

    /* End of Outputs for SubSystem: '<S104>/Convert from geodetic to  spherical coordinates ' */

    /* Memory: '<S110>/olat' */
    mainV03_56_B.olat = mainV03_56_DW.olat_PreviousInput;
  }

  /* Outputs for Enabled SubSystem: '<S104>/Convert from geodetic to  spherical coordinates ' incorporates:
   *  EnablePort: '<S108>/Enable'
   */
  if (mainV03_56_DW.Convertfromgeodetictosphericalcoordinates_MODE) {
    /* Inport: '<S108>/sp[2]' */
    mainV03_56_B.sp2 = mainV03_56_B.sincos_o1[0];

    /* Inport: '<S108>/cp[2]' */
    mainV03_56_B.cp2 = mainV03_56_B.sincos_o2[0];
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* Outputs for Iterator SubSystem: '<S108>/For Iterator Subsystem' incorporates:
       *  ForIterator: '<S155>/For Iterator'
       */
      for (s155_iter = 1; s155_iter <= mainV03_56_P.ForIterator_IterationLimit;
           s155_iter++) {
        /* Switch: '<S155>/cp[m-1] sp[m-1]' incorporates:
         *  UnitDelay: '<S155>/Unit Delay1'
         */
        if (s155_iter > mainV03_56_P.cpm1spm1_Threshold) {
          rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[0] =
            mainV03_56_DW.UnitDelay1_DSTATE[0];
          rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[1] =
            mainV03_56_DW.UnitDelay1_DSTATE[1];
        } else {
          rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[0] =
            mainV03_56_B.cp2;
          rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[1] =
            mainV03_56_B.sp2;
        }

        /* End of Switch: '<S155>/cp[m-1] sp[m-1]' */

        /* Sum: '<S155>/Sum2' incorporates:
         *  Product: '<S155>/Product1'
         *  Product: '<S155>/Product2'
         */
        rtb_IC4_idx_2 = rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[0] *
          mainV03_56_B.sp2 +
          rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[1] *
          mainV03_56_B.cp2;

        /* Assignment: '<S155>/Assignment' incorporates:
         *  Assignment: '<S155>/Assignment1'
         *  Constant: '<S155>/Constant'
         *  Constant: '<S155>/Constant1'
         */
        if (s155_iter == 1) {
          memcpy(&rtb_Assignment[0], &mainV03_56_P.Constant_Value_m[0], 11U *
                 sizeof(real_T));
          memcpy(&rtb_Assignment1[0], &mainV03_56_P.Constant1_Value_l[0], 11U *
                 sizeof(real_T));
        }

        rtb_Assignment[s155_iter - 1] = rtb_IC4_idx_2;

        /* End of Assignment: '<S155>/Assignment' */

        /* Sum: '<S155>/Sum1' incorporates:
         *  Product: '<S155>/Product3'
         *  Product: '<S155>/Product8'
         */
        rtb_Sum1_no = rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[0] *
          mainV03_56_B.cp2 -
          rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[1] *
          mainV03_56_B.sp2;

        /* Assignment: '<S155>/Assignment1' */
        rtb_Assignment1[s155_iter - 1] = rtb_Sum1_no;

        /* Update for UnitDelay: '<S155>/Unit Delay1' */
        mainV03_56_DW.UnitDelay1_DSTATE[0] = rtb_Sum1_no;
        mainV03_56_DW.UnitDelay1_DSTATE[1] = rtb_IC4_idx_2;
      }

      /* End of Outputs for SubSystem: '<S108>/For Iterator Subsystem' */
      for (s155_iter = 0; s155_iter < 11; s155_iter++) {
        /* Gain: '<S108>/Gain' */
        mainV03_56_B.Gain_k[s155_iter] = mainV03_56_P.Gain_Gain_o *
          rtb_Assignment1[s155_iter];

        /* Gain: '<S108>/Gain1' */
        mainV03_56_B.Gain1_e[s155_iter] = mainV03_56_P.Gain1_Gain_h *
          rtb_Assignment[s155_iter];
      }
    }

    /* SignalConversion: '<S108>/OutportBufferForcp[13]' incorporates:
     *  Constant: '<S108>/cp[1]'
     */
    mainV03_56_B.OutportBufferForcp13[0] = mainV03_56_P.cp1_Value;
    mainV03_56_B.OutportBufferForcp13[1] = mainV03_56_B.cp2;

    /* SignalConversion: '<S108>/OutportBufferForsp[13]' incorporates:
     *  Constant: '<S108>/sp[1]'
     */
    mainV03_56_B.OutportBufferForsp13[0] = mainV03_56_P.sp1_Value;
    mainV03_56_B.OutportBufferForsp13[1] = mainV03_56_B.sp2;

    /* SignalConversion: '<S108>/OutportBufferForcp[13]' */
    memcpy(&mainV03_56_B.OutportBufferForcp13[2], &mainV03_56_B.Gain_k[0], 11U *
           sizeof(real_T));

    /* SignalConversion: '<S108>/OutportBufferForsp[13]' */
    memcpy(&mainV03_56_B.OutportBufferForsp13[2], &mainV03_56_B.Gain1_e[0], 11U *
           sizeof(real_T));
  }

  /* End of Outputs for SubSystem: '<S104>/Convert from geodetic to  spherical coordinates ' */

  /* Saturate: '<S16>/0 to 1,000,000 m' */
  if (mainV03_56_B.plantData.LLA.Altitude_m > mainV03_56_P.uto1000000m_UpperSat)
  {
    mainV03_56_B.uto1000000m = mainV03_56_P.uto1000000m_UpperSat;
  } else if (mainV03_56_B.plantData.LLA.Altitude_m <
             mainV03_56_P.uto1000000m_LowerSat) {
    mainV03_56_B.uto1000000m = mainV03_56_P.uto1000000m_LowerSat;
  } else {
    mainV03_56_B.uto1000000m = mainV03_56_B.plantData.LLA.Altitude_m;
  }

  /* End of Saturate: '<S16>/0 to 1,000,000 m' */

  /* Gain: '<S16>/Gain' */
  mainV03_56_B.Gain = mainV03_56_P.Gain_Gain_il * mainV03_56_B.uto1000000m;
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Memory: '<S110>/oalt' */
    mainV03_56_B.oalt = mainV03_56_DW.oalt_PreviousInput;
  }

  /* Logic: '<S110>/Logical Operator' incorporates:
   *  RelationalOperator: '<S110>/Relational Operator'
   *  RelationalOperator: '<S110>/Relational Operator1'
   */
  mainV03_56_B.LogicalOperator = ((mainV03_56_B.u0deg != mainV03_56_B.olat) ||
    (mainV03_56_B.Gain != mainV03_56_B.oalt));

  /* Product: '<S109>/Product' */
  mainV03_56_B.Product_l = mainV03_56_B.sincos_o1[1] * mainV03_56_B.sincos_o1[1];

  /* Product: '<S109>/Product1' */
  mainV03_56_B.Product1_m = mainV03_56_B.sincos_o2[1] * mainV03_56_B.sincos_o2[1];
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* SignalConversion: '<S104>/HiddenBuf_InsertedFor_Convert from geodetic to  spherical coordinates_at_inport_5' */
    mainV03_56_B.HiddenBuf_InsertedFor_Convertfromgeodetictosphericalcoordinates_at_inport_5
      = mainV03_56_B.LogicalOperator;

    /* Outputs for Enabled SubSystem: '<S104>/Convert from geodetic to  spherical coordinates' incorporates:
     *  EnablePort: '<S107>/Enable'
     */
    if (rtmIsMajorTimeStep(mainV03_56_M)) {
      if (mainV03_56_B.HiddenBuf_InsertedFor_Convertfromgeodetictosphericalcoordinates_at_inport_5)
      {
        if (!mainV03_56_DW.Convertfromgeodetictosphericalcoordinates_MODE_n) {
          mainV03_56_DW.Convertfromgeodetictosphericalcoordinates_MODE_n = true;
        }
      } else {
        if (mainV03_56_DW.Convertfromgeodetictosphericalcoordinates_MODE_n) {
          mainV03_56_DW.Convertfromgeodetictosphericalcoordinates_MODE_n = false;
        }
      }
    }

    /* End of Outputs for SubSystem: '<S104>/Convert from geodetic to  spherical coordinates' */
  }

  /* Outputs for Enabled SubSystem: '<S104>/Convert from geodetic to  spherical coordinates' incorporates:
   *  EnablePort: '<S107>/Enable'
   */
  if (mainV03_56_DW.Convertfromgeodetictosphericalcoordinates_MODE_n) {
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* Product: '<S107>/b2' incorporates:
       *  Constant: '<S107>/b'
       */
      mainV03_56_B.b2 = mainV03_56_P.b_Value * mainV03_56_P.b_Value;

      /* Product: '<S107>/a2' incorporates:
       *  Constant: '<S107>/a'
       */
      mainV03_56_B.a2 = mainV03_56_P.a_Value * mainV03_56_P.a_Value;

      /* Sum: '<S150>/Sum1' */
      mainV03_56_B.c2 = mainV03_56_B.a2 - mainV03_56_B.b2;
    }

    /* Product: '<S150>/Product' */
    mainV03_56_B.Product_gg = mainV03_56_B.Product_l * mainV03_56_B.c2;

    /* Sum: '<S150>/Sum' */
    mainV03_56_B.Sum_l4 = mainV03_56_B.a2 - mainV03_56_B.Product_gg;

    /* Sqrt: '<S150>/sqrt' */
    rtb_sincos_o2_i_idx_0 = sqrt(mainV03_56_B.Sum_l4);

    /* Product: '<S107>/Product1' */
    mainV03_56_B.q1 = rtb_sincos_o2_i_idx_0 * mainV03_56_B.Gain;

    /* Product: '<S149>/Product9' */
    mainV03_56_B.Product9_i = mainV03_56_B.Product1_m * mainV03_56_B.a2;

    /* Product: '<S149>/Product10' */
    mainV03_56_B.Product10_o = mainV03_56_B.b2 * mainV03_56_B.Product_l;

    /* Sum: '<S149>/Sum7' */
    mainV03_56_B.Sum7_g = mainV03_56_B.Product9_i + mainV03_56_B.Product10_o;

    /* Sqrt: '<S149>/sqrt' */
    rtb_sincos_o2_i_idx_1 = sqrt(mainV03_56_B.Sum7_g);

    /* Sum: '<S147>/Sum8' */
    mainV03_56_B.Sum8_p = mainV03_56_B.Gain + rtb_sincos_o2_i_idx_1;

    /* Gain: '<S152>/Gain' */
    mainV03_56_B.Gain_f = mainV03_56_P.Gain_Gain_f * mainV03_56_B.q1;

    /* Product: '<S152>/Product6' */
    mainV03_56_B.Product6_p = mainV03_56_B.Gain * mainV03_56_B.Gain;
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* Product: '<S152>/a4' */
      mainV03_56_B.a4 = mainV03_56_B.a2 * mainV03_56_B.a2;

      /* Sum: '<S152>/Sum9' incorporates:
       *  Product: '<S152>/b4'
       */
      mainV03_56_B.c4 = mainV03_56_B.a4 - mainV03_56_B.b2 * mainV03_56_B.b2;
    }

    /* Product: '<S152>/Product8' */
    mainV03_56_B.Product8_d = mainV03_56_B.c4 * mainV03_56_B.Product_l;

    /* Sum: '<S152>/Sum6' */
    mainV03_56_B.Sum6_j = mainV03_56_B.a4 - mainV03_56_B.Product8_d;

    /* Product: '<S152>/Product1' */
    mainV03_56_B.Product1_oa = rtb_sincos_o2_i_idx_0 * rtb_sincos_o2_i_idx_0;

    /* Product: '<S152>/Product7' */
    mainV03_56_B.Product7_a = mainV03_56_B.Sum6_j / mainV03_56_B.Product1_oa;

    /* Sum: '<S152>/Sum5' */
    mainV03_56_B.r2 = (mainV03_56_B.Gain_f + mainV03_56_B.Product6_p) +
      mainV03_56_B.Product7_a;

    /* Sqrt: '<S107>/sqrt' */
    mainV03_56_B.sqrt_l = sqrt(mainV03_56_B.r2);

    /* Product: '<S147>/Product11' */
    mainV03_56_B.Product11_j = mainV03_56_B.Sum8_p / mainV03_56_B.sqrt_l;

    /* Sum: '<S151>/Sum2' */
    mainV03_56_B.Sum2_if = mainV03_56_B.a2 + mainV03_56_B.q1;

    /* Sum: '<S151>/Sum1' */
    mainV03_56_B.Sum1_n = mainV03_56_B.q1 + mainV03_56_B.b2;

    /* Product: '<S151>/Product1' */
    mainV03_56_B.Product1_j = mainV03_56_B.Sum1_n * mainV03_56_B.Sum1_n;

    /* Product: '<S151>/Product2' */
    mainV03_56_B.q2 = mainV03_56_B.Sum2_if * mainV03_56_B.Sum2_if /
      mainV03_56_B.Product1_j;

    /* Product: '<S148>/Product3' */
    mainV03_56_B.Product3_c = mainV03_56_B.Product1_m * mainV03_56_B.q2;

    /* Sum: '<S148>/Sum3' */
    mainV03_56_B.Sum3_c = mainV03_56_B.Product_l + mainV03_56_B.Product3_c;

    /* Product: '<S148>/Product4' incorporates:
     *  Sqrt: '<S148>/sqrt'
     */
    mainV03_56_B.Product4_l = mainV03_56_B.sincos_o1[1] / sqrt
      (mainV03_56_B.Sum3_c);

    /* Product: '<S153>/Product1' */
    mainV03_56_B.Product1_m4 = mainV03_56_B.sqrt_l * rtb_sincos_o2_i_idx_1;
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* Sum: '<S153>/Sum1' */
      mainV03_56_B.c2_e = mainV03_56_B.a2 - mainV03_56_B.b2;
    }

    /* Product: '<S153>/Product12' */
    mainV03_56_B.Product12 = mainV03_56_B.c2_e / mainV03_56_B.Product1_m4 *
      mainV03_56_B.sincos_o2[1] * mainV03_56_B.sincos_o1[1];

    /* Product: '<S154>/Product5' */
    mainV03_56_B.Product5_h = mainV03_56_B.Product4_l * mainV03_56_B.Product4_l;

    /* Sum: '<S154>/Sum4' incorporates:
     *  Constant: '<S154>/Constant'
     */
    mainV03_56_B.Sum4_g = mainV03_56_P.Constant_Value_j -
      mainV03_56_B.Product5_h;

    /* Sqrt: '<S154>/sqrt' */
    mainV03_56_B.sqrt_m = sqrt(mainV03_56_B.Sum4_g);
  }

  /* End of Outputs for SubSystem: '<S104>/Convert from geodetic to  spherical coordinates' */

  /* Product: '<S104>/aor' incorporates:
   *  Constant: '<S104>/re'
   */
  mainV03_56_B.aor = mainV03_56_P.re_Value / mainV03_56_B.sqrt_l;
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Outputs for Iterator SubSystem: '<S104>/Compute magnetic vector in spherical coordinates' incorporates:
     *  ForIterator: '<S106>/For Iterator'
     */
    /* InitializeConditions for UnitDelay: '<S106>/Unit Delay' */
    rtb_sincos_o2_i_idx_0 = mainV03_56_P.UnitDelay_InitialCondition_c;

    /* InitializeConditions for UnitDelay: '<S106>/Unit Delay2' */
    rtb_sincos_o1_n_idx_1 = mainV03_56_P.UnitDelay2_InitialCondition_j[0];
    rtb_IC4_idx_0 = mainV03_56_P.UnitDelay2_InitialCondition_j[1];
    rtb_IC4_idx_1 = mainV03_56_P.UnitDelay2_InitialCondition_j[2];
    rtb_IC3_idx_0 = mainV03_56_P.UnitDelay2_InitialCondition_j[3];
    for (s106_iter = 1; s106_iter <= mainV03_56_P.ForIterator_IterationLimit_f;
         s106_iter++) {
      /* Switch: '<S106>/ar(n)' incorporates:
       *  Product: '<S104>/ar'
       */
      if (!(s106_iter > mainV03_56_P.arn_Threshold)) {
        rtb_sincos_o2_i_idx_0 = mainV03_56_B.aor * mainV03_56_B.aor;
      }

      /* End of Switch: '<S106>/ar(n)' */

      /* Product: '<S106>/Product8' */
      rtb_sincos_o2_i_idx_0 *= mainV03_56_B.aor;

      /* Sum: '<S106>/Sum' incorporates:
       *  Constant: '<S106>/Constant'
       */
      if ((s106_iter < 0) && (mainV03_56_P.Constant_Value_io < MIN_int32_T
                              - s106_iter)) {
        rtb_PreLookUpIndexSearch_o1 = MIN_int32_T;
      } else if ((s106_iter > 0) && (mainV03_56_P.Constant_Value_io >
                  MAX_int32_T - s106_iter)) {
        rtb_PreLookUpIndexSearch_o1 = MAX_int32_T;
      } else {
        rtb_PreLookUpIndexSearch_o1 = s106_iter + mainV03_56_P.Constant_Value_io;
      }

      /* Outputs for Iterator SubSystem: '<S106>/For Iterator Subsystem' incorporates:
       *  ForIterator: '<S114>/For Iterator'
       */
      if (mainV03_56_DW.ForIterator_IterationMarker[4] != 0) {
        /* InitializeConditions for UnitDelay: '<S115>/Unit Delay1' */
        mainV03_56_DW.UnitDelay1_DSTATE_m =
          mainV03_56_P.UnitDelay1_InitialCondition_a;

        /* InitializeConditions for UnitDelay: '<S115>/Unit Delay3' */
        mainV03_56_DW.UnitDelay3_DSTATE =
          mainV03_56_P.UnitDelay3_InitialCondition;

        /* InitializeConditions for UnitDelay: '<S115>/Unit Delay2' */
        mainV03_56_DW.UnitDelay2_DSTATE_i =
          mainV03_56_P.UnitDelay2_InitialCondition;

        /* InitializeConditions for UnitDelay: '<S115>/Unit Delay4' */
        mainV03_56_DW.UnitDelay4_DSTATE =
          mainV03_56_P.UnitDelay4_InitialCondition;
      }

      for (s155_iter = 0; s155_iter < 7; s155_iter++) {
        mainV03_56_DW.ForIterator_IterationMarker[s155_iter] = 1U;
      }

      /* Sum: '<S106>/Sum' incorporates:
       *  Constant: '<S121>/Constant'
       *  Constant: '<S121>/Constant1'
       *  Constant: '<S144>/zeros(maxdef+1,maxdef+1)'
       *  Logic: '<S121>/Logical Operator'
       *  RelationalOperator: '<S121>/Relational Operator'
       *  RelationalOperator: '<S121>/Relational Operator1'
       *  SignalConversion: '<S114>/HiddenBuf_InsertedFor_Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations_at_inport_3'
       *  SignalConversion: '<S114>/HiddenBuf_InsertedFor_Time adjust the gauss coefficients_at_inport_3'
       *  Sum: '<S114>/Sum1'
       *  Switch: '<S144>/tc_old'
       *  UnitDelay: '<S144>/Unit Delay'
       */
      for (s114_iter = 1; s114_iter <= rtb_PreLookUpIndexSearch_o1; s114_iter++)
      {
        /* Sum: '<S114>/Sum1' incorporates:
         *  Constant: '<S114>/Constant'
         */
        if ((s114_iter >= 0) && (mainV03_56_P.Constant_Value_fd < s114_iter -
             MAX_int32_T)) {
          qY = MAX_int32_T;
        } else if ((s114_iter < 0) && (mainV03_56_P.Constant_Value_fd >
                    s114_iter - MIN_int32_T)) {
          qY = MIN_int32_T;
        } else {
          qY = s114_iter - mainV03_56_P.Constant_Value_fd;
        }

        /* Outputs for Enabled SubSystem: '<S114>/Time adjust the gauss coefficients' incorporates:
         *  EnablePort: '<S117>/Enable'
         */
        if (mainV03_56_B.RelationalOperator) {
          /* Outputs for Atomic SubSystem: '<S117>/If Action Subsystem' */
          /* Sum: '<S143>/Sum1' incorporates:
           *  Assignment: '<S117>/Assignment'
           *  Constant: '<S143>/Constant1'
           *  Selector: '<S143>/c[m][n]'
           *  Selector: '<S143>/cd[m][n]'
           *  Sum: '<S114>/Sum1'
           */
          if ((qY < 0) && (mainV03_56_P.Constant1_Value_e0 < MIN_int32_T - qY))
          {
            s155_iter = MIN_int32_T;
          } else if ((qY > 0) && (mainV03_56_P.Constant1_Value_e0 > MAX_int32_T
                                  - qY)) {
            s155_iter = MAX_int32_T;
          } else {
            s155_iter = qY + mainV03_56_P.Constant1_Value_e0;
          }

          rtb_Sum1_h = s155_iter - 1;

          /* End of Sum: '<S143>/Sum1' */

          /* Sum: '<S143>/Sum2' incorporates:
           *  Assignment: '<S117>/Assignment'
           *  Constant: '<S143>/Constant'
           *  Selector: '<S143>/c[m][n]'
           *  Selector: '<S143>/cd[m][n]'
           */
          if ((s106_iter < 0) && (mainV03_56_P.Constant_Value_ha < MIN_int32_T -
               s106_iter)) {
            s106_iter_0 = MIN_int32_T;
          } else if ((s106_iter > 0) && (mainV03_56_P.Constant_Value_ha >
                      MAX_int32_T - s106_iter)) {
            s106_iter_0 = MAX_int32_T;
          } else {
            s106_iter_0 = s106_iter + mainV03_56_P.Constant_Value_ha;
          }

          s155_iter = s106_iter_0 - 1;

          /* End of Sum: '<S143>/Sum2' */
          /* End of Outputs for SubSystem: '<S117>/If Action Subsystem' */

          /* Assignment: '<S117>/Assignment' incorporates:
           *  Constant: '<S117>/c[maxdef][maxdef]'
           *  Constant: '<S117>/cd[maxdef][maxdef]'
           *  Product: '<S143>/Product'
           *  Selector: '<S143>/c[m][n]'
           *  Selector: '<S143>/cd[m][n]'
           *  Sum: '<S143>/Sum'
           *  UnitDelay: '<S117>/Unit Delay'
           */
          if (mainV03_56_DW.ForIterator_IterationMarker[5] < 2) {
            mainV03_56_DW.ForIterator_IterationMarker[5] = 2U;
            memcpy(&Assignment[0], &mainV03_56_DW.UnitDelay_DSTATE_f[0], 169U *
                   sizeof(real_T));
          }

          /* Outputs for Atomic SubSystem: '<S117>/If Action Subsystem' */
          Assignment[rtb_Sum1_h + 13 * s155_iter] =
            mainV03_56_P.cdmaxdefmaxdef_Value[13 * s155_iter + rtb_Sum1_h] *
            mainV03_56_B.Sum_n + mainV03_56_P.cmaxdefmaxdef_Value[13 * s155_iter
            + rtb_Sum1_h];

          /* End of Outputs for SubSystem: '<S117>/If Action Subsystem' */
          if (s106_iter > mainV03_56_P.tc_old_Threshold) {
            memcpy(&rtb_tc_old[0], &mainV03_56_DW.UnitDelay_DSTATE_i[0], 169U *
                   sizeof(real_T));
          } else {
            memcpy(&rtb_tc_old[0], &mainV03_56_P.zerosmaxdef1maxdef1_Value[0],
                   169U * sizeof(real_T));
          }

          /* If: '<S144>/If' incorporates:
           *  Constant: '<S144>/zeros(maxdef+1,maxdef+1)'
           *  Inport: '<S146>/tc_old'
           *  Sum: '<S114>/Sum1'
           *  Switch: '<S144>/tc_old'
           *  UnitDelay: '<S144>/Unit Delay'
           */
          if (qY != 0) {
            /* Outputs for IfAction SubSystem: '<S144>/If Action Subsystem1' incorporates:
             *  ActionPort: '<S145>/Action Port'
             */
            /* Sum: '<S145>/Sum2' incorporates:
             *  Assignment: '<S145>/Assignment2'
             *  Constant: '<S145>/Constant'
             *  Selector: '<S145>/c[m][n]'
             *  Selector: '<S145>/cd[m][n]'
             */
            if ((s106_iter < 0) && (mainV03_56_P.Constant_Value_et < MIN_int32_T
                 - s106_iter)) {
              s106_iter_0 = MIN_int32_T;
            } else if ((s106_iter > 0) && (mainV03_56_P.Constant_Value_et >
                        MAX_int32_T - s106_iter)) {
              s106_iter_0 = MAX_int32_T;
            } else {
              s106_iter_0 = s106_iter + mainV03_56_P.Constant_Value_et;
            }

            rtb_Sum1_h = s106_iter_0 - 1;

            /* End of Sum: '<S145>/Sum2' */

            /* Assignment: '<S145>/Assignment2' incorporates:
             *  Constant: '<S117>/c[maxdef][maxdef]'
             *  Constant: '<S117>/cd[maxdef][maxdef]'
             *  Product: '<S145>/Product'
             *  Selector: '<S145>/c[m][n]'
             *  Selector: '<S145>/cd[m][n]'
             *  Sum: '<S145>/Sum'
             */
            if (mainV03_56_DW.ForIterator_IterationMarker[6] < 2) {
              mainV03_56_DW.ForIterator_IterationMarker[6] = 2U;
              memcpy(&Assignment2[0], &rtb_tc_old[0], 169U * sizeof(real_T));
            }

            Assignment2[rtb_Sum1_h + 13 * (qY - 1)] =
              mainV03_56_P.cdmaxdefmaxdef_Value[(qY - 1) * 13 + rtb_Sum1_h] *
              mainV03_56_B.Sum_n + mainV03_56_P.cmaxdefmaxdef_Value[(qY - 1) *
              13 + rtb_Sum1_h];

            /* Gain: '<S145>/Gain' */
            for (s155_iter = 0; s155_iter < 169; s155_iter++) {
              mainV03_56_B.Merge_l[s155_iter] = mainV03_56_P.Gain_Gain_d *
                Assignment2[s155_iter];
            }

            /* End of Gain: '<S145>/Gain' */
            /* End of Outputs for SubSystem: '<S144>/If Action Subsystem1' */
          } else {
            /* Outputs for IfAction SubSystem: '<S144>/If Action Subsystem2' incorporates:
             *  ActionPort: '<S146>/Action Port'
             */
            memcpy(&mainV03_56_B.Merge_l[0], &rtb_tc_old[0], 169U * sizeof
                   (real_T));

            /* End of Outputs for SubSystem: '<S144>/If Action Subsystem2' */
          }

          /* End of If: '<S144>/If' */

          /* Sum: '<S117>/Sum2' */
          for (s155_iter = 0; s155_iter < 169; s155_iter++) {
            mainV03_56_B.Sum2_d[s155_iter] = Assignment[s155_iter] +
              mainV03_56_B.Merge_l[s155_iter];
          }

          /* End of Sum: '<S117>/Sum2' */
        }

        /* End of Outputs for SubSystem: '<S114>/Time adjust the gauss coefficients' */

        /* Sum: '<S120>/Sum4' incorporates:
         *  Constant: '<S120>/Constant1'
         *  Constant: '<S144>/zeros(maxdef+1,maxdef+1)'
         *  SignalConversion: '<S114>/HiddenBuf_InsertedFor_Time adjust the gauss coefficients_at_inport_3'
         *  Sum: '<S114>/Sum1'
         *  Switch: '<S144>/tc_old'
         *  UnitDelay: '<S144>/Unit Delay'
         */
        rtb_Sum1_no = (real_T)qY + mainV03_56_P.Constant1_Value_p;

        /* If: '<S120>/If' incorporates:
         *  Sum: '<S114>/Sum1'
         */
        if (qY == 0) {
          /* Outputs for IfAction SubSystem: '<S120>/If Action Subsystem' incorporates:
           *  ActionPort: '<S126>/Action Port'
           */
          /* Sum: '<S126>/Sum' incorporates:
           *  Constant: '<S126>/Constant'
           */
          rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[0] = (real_T)
            s106_iter + mainV03_56_P.Constant_Value_lr;

          /* Gain: '<S126>/Gain1' incorporates:
           *  Constant: '<S126>/Constant'
           *  Product: '<S126>/Product'
           *  Selector: '<S120>/cp[m+1]'
           *  Selector: '<S126>/Selector'
           */
          rtb_sincos_o2_i_idx_1 = mainV03_56_B.Sum2_d[(((int32_T)
            rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[0] - 1) * 13 +
            (int32_T)mainV03_56_P.Constant_Value_lr) - 1] *
            mainV03_56_B.OutportBufferForcp13[(int32_T)rtb_Sum1_no - 1] *
            mainV03_56_P.Gain1_Gain;

          /* Gain: '<S126>/Gain2' incorporates:
           *  Constant: '<S126>/Constant'
           *  Product: '<S126>/Product'
           *  Selector: '<S120>/sp[m+1]'
           *  Selector: '<S126>/Selector'
           */
          rtb_IC4_idx_2 = mainV03_56_B.Sum2_d[(((int32_T)
            rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[0] - 1) * 13 +
            (int32_T)mainV03_56_P.Constant_Value_lr) - 1] *
            mainV03_56_B.OutportBufferForsp13[(int32_T)rtb_Sum1_no - 1] *
            mainV03_56_P.Gain2_Gain;

          /* End of Outputs for SubSystem: '<S120>/If Action Subsystem' */
        } else {
          /* Outputs for IfAction SubSystem: '<S120>/If Action Subsystem1' incorporates:
           *  ActionPort: '<S127>/Action Port'
           */
          /* Sum: '<S128>/Sum' incorporates:
           *  Constant: '<S128>/Constant'
           */
          rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[0] = (real_T)
            s106_iter + mainV03_56_P.Constant_Value_n;
          rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[1] = (real_T)qY +
            mainV03_56_P.Constant_Value_n;

          /* Sum: '<S127>/Sum' incorporates:
           *  Constant: '<S129>/Constant'
           *  Product: '<S127>/Product'
           *  Product: '<S127>/Product1'
           *  Selector: '<S120>/cp[m+1]'
           *  Selector: '<S120>/sp[m+1]'
           *  Selector: '<S127>/Selector'
           *  Selector: '<S127>/Selector1'
           *  Sum: '<S129>/Sum'
           */
          rtb_sincos_o2_i_idx_1 = mainV03_56_B.Sum2_d[((qY - 1) * 13 + (int32_T)
            ((real_T)s106_iter + mainV03_56_P.Constant_Value_km)) - 1] *
            mainV03_56_B.OutportBufferForsp13[(int32_T)rtb_Sum1_no - 1] +
            mainV03_56_B.Sum2_d[(((int32_T)
            rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[0] - 1) * 13 +
            (int32_T)rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[1]) -
            1] * mainV03_56_B.OutportBufferForcp13[(int32_T)rtb_Sum1_no - 1];

          /* Sum: '<S127>/Sum1' incorporates:
           *  Constant: '<S129>/Constant'
           *  Product: '<S127>/Product'
           *  Product: '<S127>/Product1'
           *  Selector: '<S120>/cp[m+1]'
           *  Selector: '<S120>/sp[m+1]'
           *  Selector: '<S127>/Selector'
           *  Selector: '<S127>/Selector1'
           *  Sum: '<S129>/Sum'
           */
          rtb_IC4_idx_2 = mainV03_56_B.Sum2_d[(((int32_T)
            rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[0] - 1) * 13 +
            (int32_T)rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[1]) -
            1] * mainV03_56_B.OutportBufferForsp13[(int32_T)rtb_Sum1_no - 1] -
            mainV03_56_B.Sum2_d[((qY - 1) * 13 + (int32_T)((real_T)s106_iter +
            mainV03_56_P.Constant_Value_km)) - 1] *
            mainV03_56_B.OutportBufferForcp13[(int32_T)rtb_Sum1_no - 1];

          /* End of Outputs for SubSystem: '<S120>/If Action Subsystem1' */
        }

        /* End of If: '<S120>/If' */

        /* Outputs for Enabled SubSystem: '<S114>/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations' incorporates:
         *  EnablePort: '<S116>/Enable'
         */
        if (mainV03_56_B.LogicalOperator) {
          /* If: '<S116>/if n == m elseif (n==1&m==0) elseif (n>1&m~=n)' incorporates:
           *  Sum: '<S114>/Sum1'
           */
          if (s106_iter == qY) {
            /* Outputs for IfAction SubSystem: '<S116>/If Action Subsystem' incorporates:
             *  ActionPort: '<S130>/Action Port'
             */
            /* Sum: '<S134>/Sum2' incorporates:
             *  Constant: '<S134>/Constant'
             */
            if ((qY >= 0) && (mainV03_56_P.Constant_Value_a < qY - MAX_int32_T))
            {
              s155_iter = MAX_int32_T;
            } else if ((qY < 0) && (mainV03_56_P.Constant_Value_a > qY -
                                    MIN_int32_T)) {
              s155_iter = MIN_int32_T;
            } else {
              s155_iter = qY - mainV03_56_P.Constant_Value_a;
            }

            /* End of Sum: '<S134>/Sum2' */

            /* Gain: '<S134>/Gain' */
            tmp = (int64_T)mainV03_56_P.Gain_Gain_ck * s155_iter;
            if (tmp > 2147483647LL) {
              tmp = 2147483647LL;
            } else {
              if (tmp < -2147483648LL) {
                tmp = -2147483648LL;
              }
            }

            rtb_Sum1_h = (int32_T)tmp;

            /* End of Gain: '<S134>/Gain' */

            /* Selector: '<S130>/Selector' incorporates:
             *  Sum: '<S134>/Sum1'
             */
            if ((s106_iter < 0) && (rtb_Sum1_h < MIN_int32_T - s106_iter)) {
              s106_iter_0 = MIN_int32_T;
            } else if ((s106_iter > 0) && (rtb_Sum1_h > MAX_int32_T - s106_iter))
            {
              s106_iter_0 = MAX_int32_T;
            } else {
              s106_iter_0 = s106_iter + rtb_Sum1_h;
            }

            /* Product: '<S130>/Product1' incorporates:
             *  Selector: '<S130>/Selector'
             *  UnitDelay: '<S116>/Unit Delay1'
             */
            mainV03_56_B.Merge_im =
              mainV03_56_DW.UnitDelay1_DSTATE_k[s106_iter_0 - 1] *
              mainV03_56_B.sqrt_m;

            /* Selector: '<S130>/Selector' incorporates:
             *  Sum: '<S134>/Sum1'
             */
            if ((s106_iter < 0) && (rtb_Sum1_h < MIN_int32_T - s106_iter)) {
              s106_iter_0 = MIN_int32_T;
            } else if ((s106_iter > 0) && (rtb_Sum1_h > MAX_int32_T - s106_iter))
            {
              s106_iter_0 = MAX_int32_T;
            } else {
              s106_iter_0 = s106_iter + rtb_Sum1_h;
            }

            /* Sum: '<S130>/Sum' incorporates:
             *  Product: '<S130>/Product'
             *  Product: '<S130>/Product2'
             *  Selector: '<S130>/Selector'
             *  Selector: '<S130>/Selector1'
             *  UnitDelay: '<S116>/Unit Delay'
             *  UnitDelay: '<S116>/Unit Delay1'
             */
            mainV03_56_B.Merge1_j = mainV03_56_DW.UnitDelay_DSTATE_k[((s106_iter
              - 1) * 13 + qY) - 1] * mainV03_56_B.sqrt_m +
              mainV03_56_DW.UnitDelay1_DSTATE_k[s106_iter_0 - 1] *
              mainV03_56_B.Product4_l;

            /* End of Outputs for SubSystem: '<S116>/If Action Subsystem' */
          } else if ((s106_iter == 1) && (qY == 0)) {
            /* Outputs for IfAction SubSystem: '<S116>/If Action Subsystem1' incorporates:
             *  ActionPort: '<S131>/Action Port'
             */
            /* Product: '<S131>/Product3' incorporates:
             *  Selector: '<S131>/Selector'
             *  UnitDelay: '<S116>/Unit Delay1'
             */
            mainV03_56_B.Merge_im = mainV03_56_DW.UnitDelay1_DSTATE_k[0] *
              mainV03_56_B.Product4_l;

            /* Sum: '<S131>/Sum' incorporates:
             *  Constant: '<S137>/Constant'
             *  Product: '<S131>/Product'
             *  Product: '<S131>/Product2'
             *  Selector: '<S131>/Selector'
             *  Selector: '<S131>/Selector1'
             *  UnitDelay: '<S116>/Unit Delay'
             *  UnitDelay: '<S116>/Unit Delay1'
             */
            mainV03_56_B.Merge1_j =
              mainV03_56_DW.UnitDelay_DSTATE_k[mainV03_56_P.Constant_Value_ci -
              1] * mainV03_56_B.Product4_l - mainV03_56_DW.UnitDelay1_DSTATE_k[0]
              * mainV03_56_B.sqrt_m;

            /* End of Outputs for SubSystem: '<S116>/If Action Subsystem1' */
          } else {
            if ((s106_iter > 1) && (s106_iter != qY)) {
              /* Outputs for IfAction SubSystem: '<S116>/If Action Subsystem2' incorporates:
               *  ActionPort: '<S132>/Action Port'
               */
              /* Sum: '<S140>/Sum' incorporates:
               *  Constant: '<S140>/Constant'
               *  Selector: '<S132>/Selector1'
               */
              if ((qY < 0) && (mainV03_56_P.Constant_Value_fe < MIN_int32_T - qY))
              {
                s155_iter = MIN_int32_T;
              } else if ((qY > 0) && (mainV03_56_P.Constant_Value_fe >
                                      MAX_int32_T - qY)) {
                s155_iter = MAX_int32_T;
              } else {
                s155_iter = qY + mainV03_56_P.Constant_Value_fe;
              }

              rtb_Sum1_h = s155_iter - 1;

              /* End of Sum: '<S140>/Sum' */

              /* Product: '<S132>/Product' incorporates:
               *  Selector: '<S132>/Selector1'
               *  UnitDelay: '<S116>/Unit Delay'
               */
              rtb_Sum1_no = mainV03_56_DW.UnitDelay_DSTATE_k[(s106_iter - 1) *
                13 + rtb_Sum1_h] * mainV03_56_B.Product4_l;

              /* Sum: '<S141>/Sum2' incorporates:
               *  Constant: '<S141>/Constant1'
               */
              if (mainV03_56_P.Constant1_Value_ge < s106_iter - MAX_int32_T) {
                rtb_Sum1_eh = MAX_int32_T;
              } else {
                rtb_Sum1_eh = s106_iter - mainV03_56_P.Constant1_Value_ge;
              }

              /* End of Sum: '<S141>/Sum2' */

              /* Switch: '<S132>/Switch' incorporates:
               *  Constant: '<S132>/Constant'
               *  DataTypeConversion: '<S141>/Data Type Conversion'
               *  RelationalOperator: '<S141>/Relational Operator'
               *  Selector: '<S132>/Selector1'
               *  UnitDelay: '<S116>/Unit Delay'
               */
              if ((rtb_Sum1_eh >= qY) > mainV03_56_P.Switch_Threshold) {
                /* Sum: '<S140>/Sum2' incorporates:
                 *  Constant: '<S140>/Constant1'
                 */
                if (mainV03_56_P.Constant1_Value_at < s106_iter - MAX_int32_T) {
                  rtb_Sum1_eh = MAX_int32_T;
                } else {
                  rtb_Sum1_eh = s106_iter - mainV03_56_P.Constant1_Value_at;
                }

                /* End of Sum: '<S140>/Sum2' */
                rtb_sincos_o1_n_idx_0 = mainV03_56_DW.UnitDelay_DSTATE_k
                  [(rtb_Sum1_eh - 1) * 13 + rtb_Sum1_h];
              } else {
                rtb_sincos_o1_n_idx_0 = mainV03_56_P.Constant_Value_i;
              }

              /* End of Switch: '<S132>/Switch' */

              /* Sum: '<S139>/Sum' incorporates:
               *  Constant: '<S139>/Constant'
               */
              if (mainV03_56_P.Constant_Value_lu > MAX_int32_T - s106_iter) {
                s155_iter = MAX_int32_T;
              } else {
                s155_iter = s106_iter + mainV03_56_P.Constant_Value_lu;
              }

              if ((qY < 0) && (mainV03_56_P.Constant_Value_lu < MIN_int32_T - qY))
              {
                s106_iter_0 = MIN_int32_T;
              } else if ((qY > 0) && (mainV03_56_P.Constant_Value_lu >
                                      MAX_int32_T - qY)) {
                s106_iter_0 = MAX_int32_T;
              } else {
                s106_iter_0 = qY + mainV03_56_P.Constant_Value_lu;
              }

              /* Gain: '<S138>/Gain' */
              tmp = (int64_T)mainV03_56_P.Gain_Gain_m * qY;
              if (tmp > 2147483647LL) {
                tmp = 2147483647LL;
              } else {
                if (tmp < -2147483648LL) {
                  tmp = -2147483648LL;
                }
              }

              rtb_Sum1_h = (int32_T)tmp;

              /* End of Gain: '<S138>/Gain' */

              /* Selector: '<S132>/Selector' incorporates:
               *  Sum: '<S138>/Sum1'
               */
              if (rtb_Sum1_h > MAX_int32_T - s106_iter) {
                rtb_Sum1_eh = MAX_int32_T;
              } else {
                rtb_Sum1_eh = s106_iter + rtb_Sum1_h;
              }

              /* Sum: '<S132>/Sum' incorporates:
               *  Constant: '<S132>/k[13][13]'
               *  Product: '<S132>/Product1'
               *  Product: '<S132>/Product4'
               *  Selector: '<S132>/Selector'
               *  Selector: '<S132>/Selector2'
               *  Sum: '<S139>/Sum'
               *  UnitDelay: '<S116>/Unit Delay1'
               */
              mainV03_56_B.Merge1_j = (rtb_Sum1_no -
                mainV03_56_DW.UnitDelay1_DSTATE_k[rtb_Sum1_eh - 1] *
                mainV03_56_B.sqrt_m) - mainV03_56_P.k1313_Value_m[((s155_iter -
                1) * 13 + s106_iter_0) - 1] * rtb_sincos_o1_n_idx_0;

              /* Sum: '<S142>/Sum2' incorporates:
               *  Constant: '<S142>/Constant1'
               */
              if (mainV03_56_P.Constant1_Value_eg < s106_iter - MAX_int32_T) {
                rtb_Sum1_eh = MAX_int32_T;
              } else {
                rtb_Sum1_eh = s106_iter - mainV03_56_P.Constant1_Value_eg;
              }

              /* End of Sum: '<S142>/Sum2' */

              /* Switch: '<S132>/Switch1' incorporates:
               *  Constant: '<S132>/Constant1'
               *  DataTypeConversion: '<S142>/Data Type Conversion'
               *  RelationalOperator: '<S142>/Relational Operator'
               *  Selector: '<S132>/Selector'
               *  UnitDelay: '<S116>/Unit Delay1'
               */
              if ((rtb_Sum1_eh >= qY) > mainV03_56_P.Switch1_Threshold) {
                /* Selector: '<S132>/Selector' incorporates:
                 *  Constant: '<S138>/Constant1'
                 *  Sum: '<S138>/Sum1'
                 *  Sum: '<S138>/Sum2'
                 */
                if (mainV03_56_P.Constant1_Value_b0 < s106_iter - MAX_int32_T) {
                  rtb_Sum1_eh = MAX_int32_T;
                } else {
                  rtb_Sum1_eh = s106_iter - mainV03_56_P.Constant1_Value_b0;
                }

                if ((rtb_Sum1_eh < 0) && (rtb_Sum1_h < MIN_int32_T - rtb_Sum1_eh))
                {
                  rtb_Sum1_eh = MIN_int32_T;
                } else if ((rtb_Sum1_eh > 0) && (rtb_Sum1_h > MAX_int32_T
                            - rtb_Sum1_eh)) {
                  rtb_Sum1_eh = MAX_int32_T;
                } else {
                  rtb_Sum1_eh += rtb_Sum1_h;
                }

                rtb_Sum1_no = mainV03_56_DW.UnitDelay1_DSTATE_k[rtb_Sum1_eh - 1];
              } else {
                rtb_Sum1_no = mainV03_56_P.Constant1_Value_d;
              }

              /* End of Switch: '<S132>/Switch1' */

              /* Selector: '<S132>/Selector' incorporates:
               *  Sum: '<S138>/Sum1'
               */
              if (rtb_Sum1_h > MAX_int32_T - s106_iter) {
                rtb_Sum1_h = MAX_int32_T;
              } else {
                rtb_Sum1_h += s106_iter;
              }

              /* Sum: '<S132>/Sum1' incorporates:
               *  Constant: '<S132>/k[13][13]'
               *  Product: '<S132>/Product2'
               *  Product: '<S132>/Product3'
               *  Selector: '<S132>/Selector'
               *  Selector: '<S132>/Selector2'
               *  Sum: '<S139>/Sum'
               *  UnitDelay: '<S116>/Unit Delay1'
               */
              mainV03_56_B.Merge_im =
                mainV03_56_DW.UnitDelay1_DSTATE_k[rtb_Sum1_h - 1] *
                mainV03_56_B.Product4_l - mainV03_56_P.k1313_Value_m[((s155_iter
                - 1) * 13 + s106_iter_0) - 1] * rtb_Sum1_no;

              /* End of Outputs for SubSystem: '<S116>/If Action Subsystem2' */
            }
          }

          /* End of If: '<S116>/if n == m elseif (n==1&m==0) elseif (n>1&m~=n)' */

          /* Sum: '<S116>/Sum' incorporates:
           *  Constant: '<S116>/Constant'
           *  Sum: '<S114>/Sum1'
           */
          if ((s106_iter < 0) && (mainV03_56_P.Constant_Value_a0 < MIN_int32_T -
               s106_iter)) {
            s155_iter = MIN_int32_T;
          } else if ((s106_iter > 0) && (mainV03_56_P.Constant_Value_a0 >
                      MAX_int32_T - s106_iter)) {
            s155_iter = MAX_int32_T;
          } else {
            s155_iter = s106_iter + mainV03_56_P.Constant_Value_a0;
          }

          if ((qY < 0) && (mainV03_56_P.Constant_Value_a0 < MIN_int32_T - qY)) {
            s106_iter_0 = MIN_int32_T;
          } else if ((qY > 0) && (mainV03_56_P.Constant_Value_a0 > MAX_int32_T -
                      qY)) {
            s106_iter_0 = MAX_int32_T;
          } else {
            s106_iter_0 = qY + mainV03_56_P.Constant_Value_a0;
          }

          /* Assignment: '<S116>/Assignment' incorporates:
           *  Sum: '<S116>/Sum'
           *  UnitDelay: '<S116>/Unit Delay'
           */
          if (mainV03_56_DW.ForIterator_IterationMarker[2] < 2) {
            mainV03_56_DW.ForIterator_IterationMarker[2] = 2U;
            memcpy(&mainV03_56_B.Assignment_i[0],
                   &mainV03_56_DW.UnitDelay_DSTATE_k[0], 169U * sizeof(real_T));
          }

          mainV03_56_B.Assignment_i[(s106_iter_0 + 13 * (s155_iter - 1)) - 1] =
            mainV03_56_B.Merge1_j;

          /* End of Assignment: '<S116>/Assignment' */

          /* Assignment: '<S116>/Assignment_snorm' incorporates:
           *  Constant: '<S133>/Constant'
           *  Gain: '<S133>/Gain'
           *  Sum: '<S116>/Sum'
           *  Sum: '<S133>/Sum1'
           *  Sum: '<S133>/Sum2'
           *  UnitDelay: '<S116>/Unit Delay1'
           */
          if (mainV03_56_DW.ForIterator_IterationMarker[3] < 2) {
            mainV03_56_DW.ForIterator_IterationMarker[3] = 2U;
            memcpy(&mainV03_56_B.Assignment_snorm[0],
                   &mainV03_56_DW.UnitDelay1_DSTATE_k[0], 169U * sizeof(real_T));
          }

          if ((s106_iter_0 >= 0) && (mainV03_56_P.Constant_Value_g0 <
               s106_iter_0 - MAX_int32_T)) {
            s106_iter_0 = MAX_int32_T;
          } else if ((s106_iter_0 < 0) && (mainV03_56_P.Constant_Value_g0 >
                      s106_iter_0 - MIN_int32_T)) {
            s106_iter_0 = MIN_int32_T;
          } else {
            s106_iter_0 -= mainV03_56_P.Constant_Value_g0;
          }

          tmp = (int64_T)mainV03_56_P.Gain_Gain_n * s106_iter_0;
          if (tmp > 2147483647LL) {
            tmp = 2147483647LL;
          } else {
            if (tmp < -2147483648LL) {
              tmp = -2147483648LL;
            }
          }

          mainV03_56_B.Assignment_snorm[(int32_T)((real_T)s155_iter + (real_T)
            (int32_T)tmp) - 1] = mainV03_56_B.Merge_im;

          /* End of Assignment: '<S116>/Assignment_snorm' */
        }

        /* End of Outputs for SubSystem: '<S114>/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations' */

        /* Selector: '<S115>/snorm[n+m*13]' incorporates:
         *  Constant: '<S119>/Constant'
         *  Gain: '<S119>/Gain'
         *  SignalConversion: '<S114>/HiddenBuf_InsertedFor_Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations_at_inport_3'
         *  Sum: '<S114>/Sum1'
         *  Sum: '<S119>/Sum1'
         */
        if ((mainV03_56_P.Constant_Value_dq < 0) && (s106_iter < MIN_int32_T
             - mainV03_56_P.Constant_Value_dq)) {
          s155_iter = MIN_int32_T;
        } else if ((mainV03_56_P.Constant_Value_dq > 0) && (s106_iter >
                    MAX_int32_T - mainV03_56_P.Constant_Value_dq)) {
          s155_iter = MAX_int32_T;
        } else {
          s155_iter = mainV03_56_P.Constant_Value_dq + s106_iter;
        }

        tmp = (int64_T)mainV03_56_P.Gain_Gain_mq * qY;
        if (tmp > 2147483647LL) {
          tmp = 2147483647LL;
        } else {
          if (tmp < -2147483648LL) {
            tmp = -2147483648LL;
          }
        }

        rtb_Sum1_h = (int32_T)tmp;
        if ((s155_iter < 0) && (rtb_Sum1_h < MIN_int32_T - s155_iter)) {
          s155_iter = MIN_int32_T;
        } else if ((s155_iter > 0) && (rtb_Sum1_h > MAX_int32_T - s155_iter)) {
          s155_iter = MAX_int32_T;
        } else {
          s155_iter += rtb_Sum1_h;
        }

        /* Product: '<S115>/par' incorporates:
         *  Selector: '<S115>/snorm[n+m*13]'
         */
        rtb_Sum1_no = mainV03_56_B.Assignment_snorm[s155_iter - 1] *
          rtb_sincos_o2_i_idx_0;

        /* Outputs for Enabled SubSystem: '<S115>/Special case - North//South Geographic Pole' incorporates:
         *  EnablePort: '<S118>/Enable'
         */
        if ((mainV03_56_B.sqrt_m == mainV03_56_P.Constant1_Value_o) &&
            (mainV03_56_P.Constant_Value_bn == qY)) {
          if (!mainV03_56_DW.SpecialcaseNorthSouthGeographicPole_MODE) {
            mainV03_56_DW.SpecialcaseNorthSouthGeographicPole_MODE = true;
          }

          /* If: '<S118>/n ==1' incorporates:
           *  Assignment: '<S123>/Assignment2'
           */
          if (s106_iter == 1) {
            /* Outputs for IfAction SubSystem: '<S118>/If Action Subsystem1' incorporates:
             *  ActionPort: '<S122>/Action Port'
             */
            /* Assignment: '<S122>/Assignment2' incorporates:
             *  Constant: '<S122>/Constant'
             *  Selector: '<S122>/pp[n-1]'
             *  Sum: '<S122>/Sum2'
             *  UnitDelay: '<S118>/Unit Delay1'
             */
            if (mainV03_56_DW.ForIterator_IterationMarker[0] < 2) {
              mainV03_56_DW.ForIterator_IterationMarker[0] = 2U;
              memcpy(&mainV03_56_B.Assignment2_e[0],
                     &mainV03_56_DW.UnitDelay1_DSTATE_n[0], 13U * sizeof(real_T));
            }

            mainV03_56_B.Assignment2_e[(int32_T)(1.0 +
              mainV03_56_P.Constant_Value_k0) - 1] =
              mainV03_56_DW.UnitDelay1_DSTATE_n[0];

            /* End of Assignment: '<S122>/Assignment2' */
            /* End of Outputs for SubSystem: '<S118>/If Action Subsystem1' */
          } else {
            /* Outputs for IfAction SubSystem: '<S118>/If Action Subsystem2' incorporates:
             *  ActionPort: '<S123>/Action Port'
             */
            if (mainV03_56_DW.ForIterator_IterationMarker[1] < 2) {
              /* Assignment: '<S123>/Assignment2' incorporates:
               *  UnitDelay: '<S118>/Unit Delay1'
               */
              mainV03_56_DW.ForIterator_IterationMarker[1] = 2U;
              memcpy(&mainV03_56_B.Assignment2_o[0],
                     &mainV03_56_DW.UnitDelay1_DSTATE_n[0], 13U * sizeof(real_T));
            }

            /* Assignment: '<S123>/Assignment2' incorporates:
             *  Constant: '<S123>/Constant'
             *  Sum: '<S123>/Sum2'
             */
            if ((s106_iter < 0) && (mainV03_56_P.Constant_Value_jj < MIN_int32_T
                 - s106_iter)) {
              s106_iter_0 = MIN_int32_T;
            } else if ((s106_iter > 0) && (mainV03_56_P.Constant_Value_jj >
                        MAX_int32_T - s106_iter)) {
              s106_iter_0 = MAX_int32_T;
            } else {
              s106_iter_0 = s106_iter + mainV03_56_P.Constant_Value_jj;
            }

            /* Sum: '<S125>/Sum' incorporates:
             *  Constant: '<S125>/Constant'
             */
            if ((qY < 0) && (mainV03_56_P.Constant_Value_ho < MIN_int32_T - qY))
            {
              s155_iter = MIN_int32_T;
            } else if ((qY > 0) && (mainV03_56_P.Constant_Value_ho > MAX_int32_T
                        - qY)) {
              s155_iter = MAX_int32_T;
            } else {
              s155_iter = qY + mainV03_56_P.Constant_Value_ho;
            }

            if ((s106_iter < 0) && (mainV03_56_P.Constant_Value_ho < MIN_int32_T
                 - s106_iter)) {
              rtb_Sum1_h = MIN_int32_T;
            } else if ((s106_iter > 0) && (mainV03_56_P.Constant_Value_ho >
                        MAX_int32_T - s106_iter)) {
              rtb_Sum1_h = MAX_int32_T;
            } else {
              rtb_Sum1_h = s106_iter + mainV03_56_P.Constant_Value_ho;
            }

            /* End of Sum: '<S125>/Sum' */

            /* Selector: '<S123>/pp[n-2] pp[n-1]' incorporates:
             *  Constant: '<S124>/Constant1'
             *  Sum: '<S124>/Sum2'
             */
            if ((s106_iter >= 0) && (mainV03_56_P.Constant1_Value_et < s106_iter
                 - MAX_int32_T)) {
              rtb_Sum1_eh = MAX_int32_T;
            } else if ((s106_iter < 0) && (mainV03_56_P.Constant1_Value_et >
                        s106_iter - MIN_int32_T)) {
              rtb_Sum1_eh = MIN_int32_T;
            } else {
              rtb_Sum1_eh = s106_iter - mainV03_56_P.Constant1_Value_et;
            }

            /* Assignment: '<S123>/Assignment2' incorporates:
             *  Constant: '<S123>/k[13][13]'
             *  Product: '<S123>/Product1'
             *  Product: '<S123>/Product2'
             *  Selector: '<S123>/Selector2'
             *  Selector: '<S123>/pp[n-2] pp[n-1]'
             *  Sum: '<S123>/Sum1'
             *  UnitDelay: '<S118>/Unit Delay1'
             */
            mainV03_56_B.Assignment2_o[s106_iter_0 - 1] =
              mainV03_56_DW.UnitDelay1_DSTATE_n[s106_iter - 1] *
              mainV03_56_B.Product4_l - mainV03_56_P.k1313_Value[((rtb_Sum1_h -
              1) * 13 + s155_iter) - 1] *
              mainV03_56_DW.UnitDelay1_DSTATE_n[rtb_Sum1_eh - 1];

            /* End of Outputs for SubSystem: '<S118>/If Action Subsystem2' */
          }

          /* End of If: '<S118>/n ==1' */

          /* SignalConversion: '<S118>/TmpSignal ConversionAtpp[n]Inport1' incorporates:
           *  UnitDelay: '<S118>/Unit Delay1'
           */
          mainV03_56_B.TmpSignalConversionAtppnInport1[0] =
            mainV03_56_DW.UnitDelay1_DSTATE_n[0];
          mainV03_56_B.TmpSignalConversionAtppnInport1[1] =
            mainV03_56_B.Assignment2_e[1];
          memcpy(&mainV03_56_B.TmpSignalConversionAtppnInport1[2],
                 &mainV03_56_B.Assignment2_o[2], 11U * sizeof(real_T));

          /* Product: '<S118>/Product2' incorporates:
           *  Constant: '<S118>/Constant'
           *  Constant: '<S118>/Constant1'
           *  Selector: '<S118>/pp[n]'
           *  Sum: '<S118>/Sum2'
           */
          mainV03_56_B.Product2_o =
            mainV03_56_B.TmpSignalConversionAtppnInport1[(int32_T)((real_T)
            s106_iter + mainV03_56_P.Constant1_Value_f) - 1] *
            rtb_sincos_o2_i_idx_0 * mainV03_56_P.Constant_Value_g *
            rtb_IC4_idx_2;
        } else {
          if (mainV03_56_DW.SpecialcaseNorthSouthGeographicPole_MODE) {
            /* Disable for Outport: '<S118>/bpp' */
            mainV03_56_B.Product2_o = mainV03_56_P.bpp_Y0;
            mainV03_56_DW.SpecialcaseNorthSouthGeographicPole_MODE = false;
          }
        }

        /* End of Outputs for SubSystem: '<S115>/Special case - North//South Geographic Pole' */

        /* Sum: '<S115>/Sum' incorporates:
         *  Constant: '<S115>/Constant'
         *  Constant: '<S121>/Constant'
         *  Constant: '<S121>/Constant1'
         *  Logic: '<S121>/Logical Operator'
         *  RelationalOperator: '<S121>/Relational Operator'
         *  RelationalOperator: '<S121>/Relational Operator1'
         *  Sum: '<S114>/Sum1'
         */
        if ((qY < 0) && (mainV03_56_P.Constant_Value_lrt < MIN_int32_T - qY)) {
          s155_iter = MIN_int32_T;
        } else if ((qY > 0) && (mainV03_56_P.Constant_Value_lrt > MAX_int32_T
                                - qY)) {
          s155_iter = MAX_int32_T;
        } else {
          s155_iter = qY + mainV03_56_P.Constant_Value_lrt;
        }

        if ((s106_iter < 0) && (mainV03_56_P.Constant_Value_lrt < MIN_int32_T
                                - s106_iter)) {
          s106_iter_0 = MIN_int32_T;
        } else if ((s106_iter > 0) && (mainV03_56_P.Constant_Value_lrt >
                    MAX_int32_T - s106_iter)) {
          s106_iter_0 = MAX_int32_T;
        } else {
          s106_iter_0 = s106_iter + mainV03_56_P.Constant_Value_lrt;
        }

        /* End of Sum: '<S115>/Sum' */

        /* Sum: '<S115>/Sum1' incorporates:
         *  Product: '<S115>/Product'
         *  Selector: '<S115>/dp[n][m]'
         *  UnitDelay: '<S115>/Unit Delay1'
         */
        mainV03_56_B.Sum1_bu = mainV03_56_DW.UnitDelay1_DSTATE_m -
          mainV03_56_B.Assignment_i[((s106_iter_0 - 1) * 13 + s155_iter) - 1] *
          rtb_sincos_o2_i_idx_1 * rtb_sincos_o2_i_idx_0;

        /* Sum: '<S115>/Sum4' incorporates:
         *  Constant: '<S115>/Constant1'
         *  Sum: '<S114>/Sum1'
         */
        if ((qY < 0) && (mainV03_56_P.Constant1_Value_dm < MIN_int32_T - qY)) {
          qY = MIN_int32_T;
        } else if ((qY > 0) && (mainV03_56_P.Constant1_Value_dm > MAX_int32_T
                                - qY)) {
          qY = MAX_int32_T;
        } else {
          qY += mainV03_56_P.Constant1_Value_dm;
        }

        /* Sum: '<S115>/Sum2' incorporates:
         *  Constant: '<S115>/fm'
         *  Product: '<S115>/Product1'
         *  Selector: '<S115>/fm[m]'
         *  UnitDelay: '<S115>/Unit Delay3'
         */
        mainV03_56_B.Sum2_p = mainV03_56_P.fm_Value[qY - 1] * rtb_Sum1_no *
          rtb_IC4_idx_2 + mainV03_56_DW.UnitDelay3_DSTATE;

        /* Sum: '<S115>/Sum4' incorporates:
         *  Constant: '<S115>/Constant1'
         */
        if ((s106_iter < 0) && (mainV03_56_P.Constant1_Value_dm < MIN_int32_T
                                - s106_iter)) {
          s106_iter_0 = MIN_int32_T;
        } else if ((s106_iter > 0) && (mainV03_56_P.Constant1_Value_dm >
                    MAX_int32_T - s106_iter)) {
          s106_iter_0 = MAX_int32_T;
        } else {
          s106_iter_0 = s106_iter + mainV03_56_P.Constant1_Value_dm;
        }

        /* Sum: '<S115>/Sum3' incorporates:
         *  Constant: '<S115>/fn'
         *  Product: '<S115>/Product2'
         *  Selector: '<S115>/fn[m]'
         *  UnitDelay: '<S115>/Unit Delay2'
         */
        mainV03_56_B.Sum3_co = mainV03_56_P.fn_Value[s106_iter_0 - 1] *
          rtb_Sum1_no * rtb_sincos_o2_i_idx_1 +
          mainV03_56_DW.UnitDelay2_DSTATE_i;

        /* Sum: '<S115>/Sum5' incorporates:
         *  UnitDelay: '<S115>/Unit Delay4'
         */
        mainV03_56_B.Sum5_k = mainV03_56_DW.UnitDelay4_DSTATE +
          mainV03_56_B.Product2_o;
        if (mainV03_56_B.RelationalOperator) {
          /* Update for Enabled SubSystem: '<S114>/Time adjust the gauss coefficients' incorporates:
           *  Update for EnablePort: '<S117>/Enable'
           */
          /* Update for UnitDelay: '<S117>/Unit Delay' */
          memcpy(&mainV03_56_DW.UnitDelay_DSTATE_f[0], &Assignment[0], 169U *
                 sizeof(real_T));

          /* Update for UnitDelay: '<S144>/Unit Delay' */
          memcpy(&mainV03_56_DW.UnitDelay_DSTATE_i[0], &mainV03_56_B.Merge_l[0],
                 169U * sizeof(real_T));

          /* End of Update for SubSystem: '<S114>/Time adjust the gauss coefficients' */
        }

        if (mainV03_56_B.LogicalOperator) {
          /* Update for Enabled SubSystem: '<S114>/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations' incorporates:
           *  Update for EnablePort: '<S116>/Enable'
           */
          /* Update for UnitDelay: '<S116>/Unit Delay' */
          memcpy(&mainV03_56_DW.UnitDelay_DSTATE_k[0],
                 &mainV03_56_B.Assignment_i[0], 169U * sizeof(real_T));

          /* Update for UnitDelay: '<S116>/Unit Delay1' */
          memcpy(&mainV03_56_DW.UnitDelay1_DSTATE_k[0],
                 &mainV03_56_B.Assignment_snorm[0], 169U * sizeof(real_T));

          /* End of Update for SubSystem: '<S114>/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations' */
        }

        /* Update for Enabled SubSystem: '<S115>/Special case - North//South Geographic Pole' incorporates:
         *  Update for EnablePort: '<S118>/Enable'
         */
        if (mainV03_56_DW.SpecialcaseNorthSouthGeographicPole_MODE) {
          /* Update for UnitDelay: '<S118>/Unit Delay1' */
          memcpy(&mainV03_56_DW.UnitDelay1_DSTATE_n[0],
                 &mainV03_56_B.TmpSignalConversionAtppnInport1[0], 13U * sizeof
                 (real_T));
        }

        /* End of Update for SubSystem: '<S115>/Special case - North//South Geographic Pole' */

        /* Update for UnitDelay: '<S115>/Unit Delay1' incorporates:
         *  SignalConversion: '<S114>/HiddenBuf_InsertedFor_Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations_at_inport_3'
         *  SignalConversion: '<S114>/HiddenBuf_InsertedFor_Time adjust the gauss coefficients_at_inport_3'
         */
        mainV03_56_DW.UnitDelay1_DSTATE_m = mainV03_56_B.Sum1_bu;

        /* Update for UnitDelay: '<S115>/Unit Delay3' */
        mainV03_56_DW.UnitDelay3_DSTATE = mainV03_56_B.Sum2_p;

        /* Update for UnitDelay: '<S115>/Unit Delay2' */
        mainV03_56_DW.UnitDelay2_DSTATE_i = mainV03_56_B.Sum3_co;

        /* Update for UnitDelay: '<S115>/Unit Delay4' */
        mainV03_56_DW.UnitDelay4_DSTATE = mainV03_56_B.Sum5_k;
      }

      /* End of Outputs for SubSystem: '<S106>/For Iterator Subsystem' */

      /* Sum: '<S106>/Sum1' incorporates:
       *  UnitDelay: '<S106>/Unit Delay2'
       */
      mainV03_56_B.Sum1_lw[0] = rtb_sincos_o1_n_idx_1 + mainV03_56_B.Sum1_bu;
      mainV03_56_B.Sum1_lw[1] = rtb_IC4_idx_0 + mainV03_56_B.Sum2_p;
      mainV03_56_B.Sum1_lw[2] = rtb_IC4_idx_1 + mainV03_56_B.Sum3_co;
      mainV03_56_B.Sum1_lw[3] = rtb_IC3_idx_0 + mainV03_56_B.Sum5_k;

      /* Update for UnitDelay: '<S106>/Unit Delay2' */
      rtb_sincos_o1_n_idx_1 = mainV03_56_B.Sum1_lw[0];
      rtb_IC4_idx_0 = mainV03_56_B.Sum1_lw[1];
      rtb_IC4_idx_1 = mainV03_56_B.Sum1_lw[2];
      rtb_IC3_idx_0 = mainV03_56_B.Sum1_lw[3];
    }

    /* End of Outputs for SubSystem: '<S104>/Compute magnetic vector in spherical coordinates' */
  }

  /* Switch: '<S158>/Switch' */
  if (mainV03_56_B.sqrt_m != 0.0) {
    /* Product: '<S158>/Product' */
    mainV03_56_B.Product_jq = mainV03_56_B.Sum1_lw[1] / mainV03_56_B.sqrt_m;
    mainV03_56_B.by = mainV03_56_B.Product_jq;
  } else {
    mainV03_56_B.by = mainV03_56_B.Sum1_lw[3];
  }

  /* End of Switch: '<S158>/Switch' */

  /* Product: '<S157>/Product1' */
  mainV03_56_B.Product1_aw = mainV03_56_B.Product11_j * mainV03_56_B.Sum1_lw[0];

  /* Product: '<S157>/Product4' */
  mainV03_56_B.Product4 = mainV03_56_B.Sum1_lw[2] * mainV03_56_B.Product12;

  /* Sum: '<S157>/Sum1' */
  mainV03_56_B.bx = (0.0 - mainV03_56_B.Product1_aw) - mainV03_56_B.Product4;

  /* Trigonometry: '<S160>/Trigonometric Function1' */
  rtb_UnaryMinus_j = rt_atan2d_snf(mainV03_56_B.by, mainV03_56_B.bx);

  /* UnitConversion: '<S162>/Unit Conversion' */
  /* Unit Conversion - from: rad to: deg
     Expression: output = (57.2958*input) + (0) */
  rtb_UnaryMinus_j *= 57.295779513082323;

  /* Product: '<S159>/Product1' */
  mainV03_56_B.Product1_h = mainV03_56_B.Product12 * mainV03_56_B.Sum1_lw[0];

  /* Product: '<S159>/Product4' */
  mainV03_56_B.Product4_m = mainV03_56_B.Sum1_lw[2] * mainV03_56_B.Product11_j;

  /* Sum: '<S159>/Sum1' */
  mainV03_56_B.bz = mainV03_56_B.Product1_h - mainV03_56_B.Product4_m;

  /* Product: '<S160>/Product' */
  mainV03_56_B.Product_k = mainV03_56_B.by * mainV03_56_B.by;

  /* Product: '<S160>/Product1' */
  mainV03_56_B.Product1_gj = mainV03_56_B.bx * mainV03_56_B.bx;

  /* Sum: '<S160>/Sum' */
  mainV03_56_B.Sum_pr = mainV03_56_B.Product_k + mainV03_56_B.Product1_gj;

  /* Sqrt: '<S160>/sqrt1' */
  rtb_UnaryMinus_o = sqrt(mainV03_56_B.Sum_pr);

  /* Trigonometry: '<S160>/Trigonometric Function' */
  rtb_UnaryMinus_o = rt_atan2d_snf(mainV03_56_B.bz, rtb_UnaryMinus_o);

  /* UnitConversion: '<S161>/Unit Conversion' */
  /* Unit Conversion - from: rad to: deg
     Expression: output = (57.2958*input) + (0) */
  rtb_UnaryMinus_o *= 57.295779513082323;

  /* UnitConversion: '<S105>/Unit Conversion' */
  /* Unit Conversion - from: deg to: rad
     Expression: output = (0.0174533*input) + (0) */
  rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[0] = 0.017453292519943295
    * rtb_UnaryMinus_j;
  rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[1] = 0.017453292519943295
    * rtb_UnaryMinus_o;

  /* Trigonometry: '<S99>/sincos' */
  rtb_Sum1_no = cos(rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[0]);
  rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[0] = sin
    (rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[0]);

  /* Product: '<S160>/Product2' */
  mainV03_56_B.Product2_n = mainV03_56_B.bz * mainV03_56_B.bz;

  /* Sum: '<S160>/Sum1' */
  mainV03_56_B.Sum1_p = mainV03_56_B.Product2_n + mainV03_56_B.Sum_pr;

  /* Sqrt: '<S160>/sqrt' */
  rtb_UnaryMinus_j = sqrt(mainV03_56_B.Sum1_p);

  /* Product: '<S99>/h1' incorporates:
   *  Trigonometry: '<S99>/sincos'
   */
  mainV03_56_B.h1 = cos(rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[1])
    * rtb_UnaryMinus_j;

  /* Product: '<S99>/x1' */
  mainV03_56_B.x1 = rtb_Sum1_no * mainV03_56_B.h1;

  /* Product: '<S99>/y1' */
  mainV03_56_B.y1 = rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[0] *
    mainV03_56_B.h1;

  /* Product: '<S99>/z1' incorporates:
   *  Trigonometry: '<S99>/sincos'
   */
  mainV03_56_B.z1 = sin(rtb_TmpSignalConversionAtPreLookUpIndexSearchInport2[1])
    * rtb_UnaryMinus_j;

  /* SignalConversion: '<S2>/TmpSignal ConversionAtenvDataBusInport12' */
  mainV03_56_B.envData.MagneticField_ned[0] = mainV03_56_B.x1;
  mainV03_56_B.envData.MagneticField_ned[1] = mainV03_56_B.y1;
  mainV03_56_B.envData.MagneticField_ned[2] = mainV03_56_B.z1;

  /* Product: '<S2>/To Body Axes' */
  for (rtb_Sum1_eh = 0; rtb_Sum1_eh < 3; rtb_Sum1_eh++) {
    mainV03_56_B.envData.Fgravity[rtb_Sum1_eh] = 0.0;
    mainV03_56_B.envData.Fgravity[rtb_Sum1_eh] +=
      mainV03_56_B.plantData.DCM_body_earth[rtb_Sum1_eh] *
      mainV03_56_B.MathFunction[0];
    mainV03_56_B.envData.Fgravity[rtb_Sum1_eh] +=
      mainV03_56_B.plantData.DCM_body_earth[rtb_Sum1_eh + 3] *
      mainV03_56_B.MathFunction[1];
    mainV03_56_B.envData.Fgravity[rtb_Sum1_eh] +=
      mainV03_56_B.plantData.DCM_body_earth[rtb_Sum1_eh + 6] *
      mainV03_56_B.MathFunction[2];
  }

  /* End of Product: '<S2>/To Body Axes' */

  /* Clock: '<S47>/Time' */
  rtb_UnaryMinus_j = mainV03_56_M->Timing.t[0];

  /* RelationalOperator: '<S47>/Gust Start' incorporates:
   *  Constant: '<S47>/Gust start time'
   */
  rtb_UnaryMinus_j = (rtb_UnaryMinus_j >= mainV03_56_P.DiscreteWindGustModel_t_0);

  /* Logic: '<S47>/Logical Operator2' incorporates:
   *  Constant: '<S47>/Constant'
   */
  mainV03_56_B.LogicalOperator2 = ((rtb_UnaryMinus_j != 0.0) &&
    (mainV03_56_P.DiscreteWindGustModel_Gx != 0.0));
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* SignalConversion: '<S47>/HiddenBuf_InsertedFor_Distance into gust (x)_at_inport_1' */
    mainV03_56_B.HiddenBuf_InsertedFor_Distanceintogustx_at_inport_1 =
      mainV03_56_B.LogicalOperator2;

    /* Outputs for Enabled SubSystem: '<S47>/Distance into gust (x)' incorporates:
     *  EnablePort: '<S50>/Enable'
     */
    if (rtmIsMajorTimeStep(mainV03_56_M)) {
      if (mainV03_56_B.HiddenBuf_InsertedFor_Distanceintogustx_at_inport_1) {
        if (!mainV03_56_DW.Distanceintogustx_MODE) {
          /* InitializeConditions for Integrator: '<S50>/Distance into Gust (x) (Limited to gust length d)' */
          mainV03_56_X.DistanceintoGustxLimitedtogustlengthd_CSTATE_a =
            mainV03_56_P.DistanceintoGustxLimitedtogustlengthd_IC;
          mainV03_56_DW.Distanceintogustx_MODE = true;
        }
      } else {
        if (mainV03_56_DW.Distanceintogustx_MODE) {
          mainV03_56_DW.Distanceintogustx_MODE = false;
        }
      }
    }

    /* End of Outputs for SubSystem: '<S47>/Distance into gust (x)' */
  }

  /* Outputs for Enabled SubSystem: '<S47>/Distance into gust (x)' incorporates:
   *  EnablePort: '<S50>/Enable'
   */
  if (mainV03_56_DW.Distanceintogustx_MODE) {
    /* Integrator: '<S50>/Distance into Gust (x) (Limited to gust length d)' */
    /* Limited  Integrator  */
    if (mainV03_56_X.DistanceintoGustxLimitedtogustlengthd_CSTATE_a >=
        mainV03_56_P.Distanceintogustx_d_m) {
      mainV03_56_X.DistanceintoGustxLimitedtogustlengthd_CSTATE_a =
        mainV03_56_P.Distanceintogustx_d_m;
    } else {
      if (mainV03_56_X.DistanceintoGustxLimitedtogustlengthd_CSTATE_a <=
          mainV03_56_P.DistanceintoGustxLimitedtogustlengthd_LowerSat) {
        mainV03_56_X.DistanceintoGustxLimitedtogustlengthd_CSTATE_a =
          mainV03_56_P.DistanceintoGustxLimitedtogustlengthd_LowerSat;
      }
    }

    mainV03_56_B.DistanceintoGustxLimitedtogustlengthd =
      mainV03_56_X.DistanceintoGustxLimitedtogustlengthd_CSTATE_a;

    /* End of Integrator: '<S50>/Distance into Gust (x) (Limited to gust length d)' */
  }

  /* End of Outputs for SubSystem: '<S47>/Distance into gust (x)' */

  /* Logic: '<S47>/Logical Operator1' incorporates:
   *  Constant: '<S47>/Constant1'
   */
  mainV03_56_B.LogicalOperator1 = ((rtb_UnaryMinus_j != 0.0) &&
    (mainV03_56_P.DiscreteWindGustModel_Gy != 0.0));
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* SignalConversion: '<S47>/HiddenBuf_InsertedFor_Distance into gust (y)_at_inport_1' */
    mainV03_56_B.HiddenBuf_InsertedFor_Distanceintogusty_at_inport_1 =
      mainV03_56_B.LogicalOperator1;
  }

  /* Outputs for Enabled SubSystem: '<S47>/Distance into gust (y)' */
  mainV03_56_Distanceintogusty(mainV03_56_M,
    mainV03_56_B.HiddenBuf_InsertedFor_Distanceintogusty_at_inport_1,
    &mainV03_56_B.Distanceintogusty, &mainV03_56_DW.Distanceintogusty,
    (P_Distanceintogusty_mainV03_56_T *)&mainV03_56_P.Distanceintogusty,
    &mainV03_56_X.Distanceintogusty, mainV03_56_P.Distanceintogusty_d_m);

  /* End of Outputs for SubSystem: '<S47>/Distance into gust (y)' */

  /* Logic: '<S47>/Logical Operator3' incorporates:
   *  Constant: '<S47>/Constant2'
   */
  mainV03_56_B.LogicalOperator3 = ((rtb_UnaryMinus_j != 0.0) &&
    (mainV03_56_P.DiscreteWindGustModel_Gz != 0.0));
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* SignalConversion: '<S47>/HiddenBuf_InsertedFor_Distance into gust (z)_at_inport_1' */
    mainV03_56_B.HiddenBuf_InsertedFor_Distanceintogustz_at_inport_1 =
      mainV03_56_B.LogicalOperator3;
  }

  /* Outputs for Enabled SubSystem: '<S47>/Distance into gust (z)' */
  mainV03_56_Distanceintogusty(mainV03_56_M,
    mainV03_56_B.HiddenBuf_InsertedFor_Distanceintogustz_at_inport_1,
    &mainV03_56_B.Distanceintogustz, &mainV03_56_DW.Distanceintogustz,
    (P_Distanceintogusty_mainV03_56_T *)&mainV03_56_P.Distanceintogustz,
    &mainV03_56_X.Distanceintogustz, mainV03_56_P.Distanceintogustz_d_m);

  /* End of Outputs for SubSystem: '<S47>/Distance into gust (z)' */

  /* Gain: '<S47>/pi//Gust length' */
  mainV03_56_B.piGustlength[0] = 3.1415926535897931 /
    mainV03_56_P.DiscreteWindGustModel_d_m[0] *
    mainV03_56_B.DistanceintoGustxLimitedtogustlengthd;
  mainV03_56_B.piGustlength[1] = 3.1415926535897931 /
    mainV03_56_P.DiscreteWindGustModel_d_m[1] *
    mainV03_56_B.Distanceintogusty.DistanceintoGustxLimitedtogustlengthd;
  mainV03_56_B.piGustlength[2] = 3.1415926535897931 /
    mainV03_56_P.DiscreteWindGustModel_d_m[2] *
    mainV03_56_B.Distanceintogustz.DistanceintoGustxLimitedtogustlengthd;

  /* Sum: '<S47>/Sum1' incorporates:
   *  Constant: '<S47>/2'
   *  Trigonometry: '<S47>/cos(pi*x//dm)'
   */
  mainV03_56_B.Sum1_m[0] = mainV03_56_P.u_Value - cos(mainV03_56_B.piGustlength
    [0]);

  /* Gain: '<S47>/Gust magnitude//2.0' */
  mainV03_56_B.Gustmagnitude20[0] = mainV03_56_P.DiscreteWindGustModel_v_m[0] /
    2.0 * mainV03_56_B.Sum1_m[0];

  /* Sum: '<S47>/Sum1' incorporates:
   *  Constant: '<S47>/2'
   *  Trigonometry: '<S47>/cos(pi*x//dm)'
   */
  mainV03_56_B.Sum1_m[1] = mainV03_56_P.u_Value - cos(mainV03_56_B.piGustlength
    [1]);

  /* Gain: '<S47>/Gust magnitude//2.0' */
  mainV03_56_B.Gustmagnitude20[1] = mainV03_56_P.DiscreteWindGustModel_v_m[1] /
    2.0 * mainV03_56_B.Sum1_m[1];

  /* Sum: '<S47>/Sum1' incorporates:
   *  Constant: '<S47>/2'
   *  Trigonometry: '<S47>/cos(pi*x//dm)'
   */
  mainV03_56_B.Sum1_m[2] = mainV03_56_P.u_Value - cos(mainV03_56_B.piGustlength
    [2]);

  /* Gain: '<S47>/Gust magnitude//2.0' */
  mainV03_56_B.Gustmagnitude20[2] = mainV03_56_P.DiscreteWindGustModel_v_m[2] /
    2.0 * mainV03_56_B.Sum1_m[2];

  /* DotProduct: '<S15>/Dot Product' */
  rtb_UnaryMinus_j = (mainV03_56_B.plantData.V_ned[0] *
                      mainV03_56_B.plantData.V_ned[0] +
                      mainV03_56_B.plantData.V_ned[1] *
                      mainV03_56_B.plantData.V_ned[1]) +
    mainV03_56_B.plantData.V_ned[2] * mainV03_56_B.plantData.V_ned[2];
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* UnitConversion: '<S54>/Unit Conversion' incorporates:
     *  Constant: '<S48>/Wind direction'
     */
    /* Unit Conversion - from: deg to: rad
       Expression: output = (0.0174533*input) + (0) */
    mainV03_56_B.UnitConversion = 0.017453292519943295 *
      mainV03_56_P.DrydenWindTurbulenceModelDiscreteqr_Wdeg;
  }

  /* UnitConversion: '<S57>/Unit Conversion' */
  /* Unit Conversion - from: m to: ft
     Expression: output = (3.28084*input) + (0) */
  rtb_sincos_o2_i_idx_0 = 3.280839895013123 *
    mainV03_56_B.plantData.LLA.Altitude_m;

  /* Saturate: '<S91>/Limit Function 10ft to 1000ft' */
  if (rtb_sincos_o2_i_idx_0 > mainV03_56_P.LimitFunction10ftto1000ft_UpperSat) {
    mainV03_56_B.LimitFunction10ftto1000ft =
      mainV03_56_P.LimitFunction10ftto1000ft_UpperSat;
  } else if (rtb_sincos_o2_i_idx_0 <
             mainV03_56_P.LimitFunction10ftto1000ft_LowerSat) {
    mainV03_56_B.LimitFunction10ftto1000ft =
      mainV03_56_P.LimitFunction10ftto1000ft_LowerSat;
  } else {
    mainV03_56_B.LimitFunction10ftto1000ft = rtb_sincos_o2_i_idx_0;
  }

  /* End of Saturate: '<S91>/Limit Function 10ft to 1000ft' */
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* UnitConversion: '<S93>/Unit Conversion' incorporates:
     *  Constant: '<S92>/Medium//High Altitude'
     */
    /* Unit Conversion - from: m to: ft
       Expression: output = (3.28084*input) + (0) */
    mainV03_56_B.UnitConversion_j = 3.280839895013123 *
      mainV03_56_P.DrydenWindTurbulenceModelDiscreteqr_L_high;

    /* UnitConversion: '<S64>/Unit Conversion' incorporates:
     *  Constant: '<S48>/Windspeed at 20ft (6m)'
     */
    /* Unit Conversion - from: m/s to: ft/s
       Expression: output = (3.28084*input) + (0) */
    rtb_Switch = 3.280839895013123 *
      mainV03_56_P.DrydenWindTurbulenceModelDiscreteqr_W20;

    /* Gain: '<S74>/sigma_wg ' */
    mainV03_56_B.sigma_wg = mainV03_56_P.sigma_wg_Gain * rtb_Switch;
  }

  /* Gain: '<S62>/Lw' */
  mainV03_56_B.Lw[0] = mainV03_56_P.Lw_Gain *
    mainV03_56_B.LimitFunction10ftto1000ft;
  mainV03_56_B.Lw[1] = mainV03_56_P.Lw_Gain * mainV03_56_B.UnitConversion_j;

  /* PreLookup: '<S73>/PreLook-Up Index Search  (altitude)' */
  mainV03_56_B.PreLookUpIndexSearchaltitude_o1 = plook_bincpa
    (rtb_sincos_o2_i_idx_0,
     mainV03_56_P.PreLookUpIndexSearchaltitude_BreakpointsData, 11U,
     &mainV03_56_B.PreLookUpIndexSearchaltitude_o2,
     &mainV03_56_DW.PreLookUpIndexSearchaltitude_DWORK1);
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* PreLookup: '<S73>/PreLook-Up Index Search  (prob of exceed)' incorporates:
     *  Constant: '<S73>/Probability of  Exceedance'
     */
    mainV03_56_B.PreLookUpIndexSearchprobofexceed_o1 = plook_bincpa
      (mainV03_56_P.DrydenWindTurbulenceModelDiscreteqr_TurbProb,
       mainV03_56_P.PreLookUpIndexSearchprobofexceed_BreakpointsData, 6U,
       &mainV03_56_B.PreLookUpIndexSearchprobofexceed_o2,
       &mainV03_56_DW.PreLookUpIndexSearchprobofexceed_DWORK1);
  }

  /* Interpolation_n-D: '<S73>/Medium//High Altitude Intensity' */
  frac[0] = mainV03_56_B.PreLookUpIndexSearchaltitude_o2;
  frac[1] = mainV03_56_B.PreLookUpIndexSearchprobofexceed_o2;
  bpIndex[0] = mainV03_56_B.PreLookUpIndexSearchaltitude_o1;
  bpIndex[1] = mainV03_56_B.PreLookUpIndexSearchprobofexceed_o1;
  mainV03_56_B.MediumHighAltitudeIntensity = intrp2d_la_pw(bpIndex, frac,
    mainV03_56_P.MediumHighAltitudeIntensity_Table, 12U,
    mainV03_56_P.MediumHighAltitudeIntensity_maxIndex);
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[2] == 0) {
    /* Gain: '<S65>/Output' incorporates:
     *  RandomNumber: '<S65>/White Noise'
     */
    mainV03_56_B.Output = mainV03_56_P.Output_Gain * mainV03_56_DW.NextOutput;
  }

  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* UnitConversion: '<S58>/Unit Conversion' incorporates:
     *  Constant: '<S48>/Wingspan'
     */
    /* Unit Conversion - from: m to: ft
       Expression: output = (3.28084*input) + (0) */
    mainV03_56_B.UnitConversion_l = 3.280839895013123 *
      mainV03_56_P.DrydenWindTurbulenceModelDiscreteqr_Wingspan;

    /* Outputs for Enabled SubSystem: '<S55>/Hpgw' incorporates:
     *  EnablePort: '<S67>/Enable'
     */
    if (rtmIsMajorTimeStep(mainV03_56_M)) {
      /* Constant: '<S55>/Constant1' */
      if (mainV03_56_P.DrydenWindTurbulenceModelDiscreteqr_T_on > 0.0) {
        if (!mainV03_56_DW.Hpgw_MODE) {
          /* InitializeConditions for UnitDelay: '<S67>/Unit Delay' */
          mainV03_56_DW.UnitDelay_DSTATE_c[0] =
            mainV03_56_P.UnitDelay_InitialCondition;
          mainV03_56_DW.UnitDelay_DSTATE_c[1] =
            mainV03_56_P.UnitDelay_InitialCondition;
          mainV03_56_DW.Hpgw_MODE = true;
        }
      } else {
        if (mainV03_56_DW.Hpgw_MODE) {
          /* Disable for Outport: '<S67>/pgw' */
          mainV03_56_B.Sum_ax[0] = mainV03_56_P.pgw_Y0;
          mainV03_56_B.Sum_ax[1] = mainV03_56_P.pgw_Y0;
          mainV03_56_DW.Hpgw_MODE = false;
        }
      }

      /* End of Constant: '<S55>/Constant1' */
    }

    /* End of Outputs for SubSystem: '<S55>/Hpgw' */
  }

  /* Outputs for Enabled SubSystem: '<S55>/Hpgw' incorporates:
   *  EnablePort: '<S67>/Enable'
   */
  if (mainV03_56_DW.Hpgw_MODE) {
    /* Product: '<S67>/w2' */
    mainV03_56_B.w2_bi[0] = mainV03_56_B.Lw[0] * mainV03_56_B.UnitConversion_l;

    /* Product: '<S67>/w1' incorporates:
     *  Constant: '<S67>/Constant2'
     *  Sqrt: '<S67>/sqrt'
     */
    mainV03_56_B.ap[0] = mainV03_56_P.Constant2_Value / sqrt(mainV03_56_B.w2_bi
      [0]);

    /* Gain: '<S67>/2' */
    mainV03_56_B.u_j[0] = mainV03_56_P.u_Gain * mainV03_56_B.ap[0];

    /* Product: '<S67>/w4' */
    mainV03_56_B.w4[0] = mainV03_56_B.w2_bi[0] * mainV03_56_B.UnitConversion_l;

    /* Product: '<S67>/w2' */
    mainV03_56_B.w2_bi[1] = mainV03_56_B.Lw[1] * mainV03_56_B.UnitConversion_l;

    /* Product: '<S67>/w1' incorporates:
     *  Constant: '<S67>/Constant2'
     *  Sqrt: '<S67>/sqrt'
     */
    mainV03_56_B.ap[1] = mainV03_56_P.Constant2_Value / sqrt(mainV03_56_B.w2_bi
      [1]);

    /* Gain: '<S67>/2' */
    mainV03_56_B.u_j[1] = mainV03_56_P.u_Gain * mainV03_56_B.ap[1];

    /* Product: '<S67>/w4' */
    mainV03_56_B.w4[1] = mainV03_56_B.w2_bi[1] * mainV03_56_B.UnitConversion_l;

    /* Math: '<S67>/Math Function' incorporates:
     *  Constant: '<S67>/Constant3'
     */
    if ((mainV03_56_B.w4[0] < 0.0) && (mainV03_56_P.Constant3_Value > floor
         (mainV03_56_P.Constant3_Value))) {
      riseValLimit = -rt_powd_snf(-mainV03_56_B.w4[0],
        mainV03_56_P.Constant3_Value);
    } else {
      riseValLimit = rt_powd_snf(mainV03_56_B.w4[0],
        mainV03_56_P.Constant3_Value);
    }

    /* Product: '<S67>/w3' incorporates:
     *  Constant: '<S67>/Constant1'
     */
    mainV03_56_B.sp[0] = mainV03_56_P.Constant1_Value / riseValLimit *
      mainV03_56_B.sigma_wg;

    /* Math: '<S67>/Math Function' incorporates:
     *  Constant: '<S67>/Constant3'
     */
    if ((mainV03_56_B.w4[1] < 0.0) && (mainV03_56_P.Constant3_Value > floor
         (mainV03_56_P.Constant3_Value))) {
      riseValLimit = -rt_powd_snf(-mainV03_56_B.w4[1],
        mainV03_56_P.Constant3_Value);
    } else {
      riseValLimit = rt_powd_snf(mainV03_56_B.w4[1],
        mainV03_56_P.Constant3_Value);
    }

    /* Product: '<S67>/w3' incorporates:
     *  Constant: '<S67>/Constant1'
     */
    mainV03_56_B.sp[1] = mainV03_56_P.Constant1_Value / riseValLimit *
      mainV03_56_B.MediumHighAltitudeIntensity;

    /* Product: '<S67>/Lug//V1' incorporates:
     *  Sqrt: '<S67>/sqrt1'
     */
    mainV03_56_B.LugV1_f[0] = sqrt(mainV03_56_B.u_j[0]) * mainV03_56_B.sp[0] *
      mainV03_56_B.Output;

    /* Gain: '<S67>/dt' */
    mainV03_56_B.dt_d[0] = mainV03_56_P.dt_Gain * mainV03_56_B.ap[0];

    /* Sum: '<S67>/Sum1' incorporates:
     *  Constant: '<S67>/Constant'
     */
    mainV03_56_B.Sum1_cw[0] = mainV03_56_P.Constant_Value - mainV03_56_B.dt_d[0];

    /* Product: '<S67>/Lug//V1' incorporates:
     *  Sqrt: '<S67>/sqrt1'
     */
    mainV03_56_B.LugV1_f[1] = sqrt(mainV03_56_B.u_j[1]) * mainV03_56_B.sp[1] *
      mainV03_56_B.Output;

    /* Gain: '<S67>/dt' */
    mainV03_56_B.dt_d[1] = mainV03_56_P.dt_Gain * mainV03_56_B.ap[1];

    /* Sum: '<S67>/Sum1' incorporates:
     *  Constant: '<S67>/Constant'
     */
    mainV03_56_B.Sum1_cw[1] = mainV03_56_P.Constant_Value - mainV03_56_B.dt_d[1];
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[2] == 0) {
      /* UnitDelay: '<S67>/Unit Delay' */
      mainV03_56_B.UnitDelay_o[0] = mainV03_56_DW.UnitDelay_DSTATE_c[0];
      mainV03_56_B.UnitDelay_o[1] = mainV03_56_DW.UnitDelay_DSTATE_c[1];
    }

    /* Product: '<S67>/Lug//V2' */
    mainV03_56_B.LugV2_pf[0] = mainV03_56_B.Sum1_cw[0] *
      mainV03_56_B.UnitDelay_o[0];

    /* Sum: '<S67>/Sum' */
    mainV03_56_B.Sum_ax[0] = mainV03_56_B.LugV2_pf[0] + mainV03_56_B.LugV1_f[0];

    /* Product: '<S67>/Lug//V2' */
    mainV03_56_B.LugV2_pf[1] = mainV03_56_B.Sum1_cw[1] *
      mainV03_56_B.UnitDelay_o[1];

    /* Sum: '<S67>/Sum' */
    mainV03_56_B.Sum_ax[1] = mainV03_56_B.LugV2_pf[1] + mainV03_56_B.LugV1_f[1];
  }

  /* End of Outputs for SubSystem: '<S55>/Hpgw' */

  /* Sqrt: '<S15>/vt' */
  mainV03_56_B.Airspeed = sqrt(rtb_UnaryMinus_j);

  /* UnitConversion: '<S63>/Unit Conversion' */
  /* Unit Conversion - from: m/s to: ft/s
     Expression: output = (3.28084*input) + (0) */
  rtb_sincos_o2_i_idx_1 = 3.280839895013123 * mainV03_56_B.Airspeed;
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[2] == 0) {
    /* Gain: '<S66>/Output' incorporates:
     *  RandomNumber: '<S66>/White Noise'
     */
    mainV03_56_B.Output_o[0] = mainV03_56_P.Output_Gain_g[0] *
      mainV03_56_DW.NextOutput_a[0];
    mainV03_56_B.Output_o[1] = mainV03_56_P.Output_Gain_g[1] *
      mainV03_56_DW.NextOutput_a[1];
    mainV03_56_B.Output_o[2] = mainV03_56_P.Output_Gain_g[2] *
      mainV03_56_DW.NextOutput_a[2];
  }

  /* Outputs for Enabled SubSystem: '<S55>/Hqgw' incorporates:
   *  EnablePort: '<S68>/Enable'
   */
  /* Outputs for Enabled SubSystem: '<S56>/Hwgw(z)' incorporates:
   *  EnablePort: '<S72>/Enable'
   */
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    if (rtmIsMajorTimeStep(mainV03_56_M)) {
      /* Constant: '<S56>/Constant3' */
      if (mainV03_56_P.DrydenWindTurbulenceModelDiscreteqr_T_on > 0.0) {
        if (!mainV03_56_DW.Hwgwz_MODE) {
          /* InitializeConditions for UnitDelay: '<S72>/Unit Delay' */
          mainV03_56_DW.UnitDelay_DSTATE_l[0] =
            mainV03_56_P.UnitDelay_InitialCondition_kv;
          mainV03_56_DW.UnitDelay_DSTATE_l[1] =
            mainV03_56_P.UnitDelay_InitialCondition_kv;
          mainV03_56_DW.Hwgwz_MODE = true;
        }
      } else {
        if (mainV03_56_DW.Hwgwz_MODE) {
          /* Disable for Outport: '<S72>/wgw' */
          mainV03_56_B.Sum_da[0] = mainV03_56_P.wgw_Y0;
          mainV03_56_B.Sum_da[1] = mainV03_56_P.wgw_Y0;
          mainV03_56_DW.Hwgwz_MODE = false;
        }
      }
    }

    if (rtmIsMajorTimeStep(mainV03_56_M)) {
      /* Constant: '<S55>/Constant2' */
      if (mainV03_56_P.DrydenWindTurbulenceModelDiscreteqr_T_on > 0.0) {
        if (!mainV03_56_DW.Hqgw_MODE) {
          /* InitializeConditions for UnitDelay: '<S68>/Unit Delay' */
          mainV03_56_DW.UnitDelay_DSTATE_fj[0] =
            mainV03_56_P.UnitDelay_InitialCondition_f;

          /* InitializeConditions for UnitDelay: '<S68>/Unit Delay1' */
          mainV03_56_DW.UnitDelay1_DSTATE_j[0] =
            mainV03_56_P.UnitDelay1_InitialCondition;

          /* InitializeConditions for UnitDelay: '<S68>/Unit Delay' */
          mainV03_56_DW.UnitDelay_DSTATE_fj[1] =
            mainV03_56_P.UnitDelay_InitialCondition_f;

          /* InitializeConditions for UnitDelay: '<S68>/Unit Delay1' */
          mainV03_56_DW.UnitDelay1_DSTATE_j[1] =
            mainV03_56_P.UnitDelay1_InitialCondition;
          mainV03_56_DW.Hqgw_MODE = true;
        }
      } else {
        if (mainV03_56_DW.Hqgw_MODE) {
          /* Disable for Outport: '<S68>/qgw' */
          mainV03_56_B.Sum1_gx[0] = mainV03_56_P.qgw_Y0;
          mainV03_56_B.Sum1_gx[1] = mainV03_56_P.qgw_Y0;
          mainV03_56_DW.Hqgw_MODE = false;
        }
      }

      /* End of Constant: '<S55>/Constant2' */
    }
  }

  /* End of Outputs for SubSystem: '<S55>/Hqgw' */
  if (mainV03_56_DW.Hwgwz_MODE) {
    /* Product: '<S72>/V//Lwg' */
    mainV03_56_B.VLwg[0] = rtb_sincos_o2_i_idx_1 / mainV03_56_B.Lw[0];

    /* Gain: '<S72>/2' */
    mainV03_56_B.u[0] = mainV03_56_P.u_Gain_h * mainV03_56_B.VLwg[0];

    /* Gain: '<S72>/dt' */
    mainV03_56_B.dt[0] = mainV03_56_P.dt_Gain_o * mainV03_56_B.VLwg[0];

    /* Sum: '<S72>/Sum1' incorporates:
     *  Constant: '<S72>/Constant'
     */
    mainV03_56_B.Sum1_a[0] = mainV03_56_P.Constant_Value_o - mainV03_56_B.dt[0];

    /* Product: '<S72>/V//Lwg' */
    mainV03_56_B.VLwg[1] = rtb_sincos_o2_i_idx_1 / mainV03_56_B.Lw[1];

    /* Gain: '<S72>/2' */
    mainV03_56_B.u[1] = mainV03_56_P.u_Gain_h * mainV03_56_B.VLwg[1];

    /* Gain: '<S72>/dt' */
    mainV03_56_B.dt[1] = mainV03_56_P.dt_Gain_o * mainV03_56_B.VLwg[1];

    /* Sum: '<S72>/Sum1' incorporates:
     *  Constant: '<S72>/Constant'
     */
    mainV03_56_B.Sum1_a[1] = mainV03_56_P.Constant_Value_o - mainV03_56_B.dt[1];

    /* Product: '<S72>/Lug//V1' incorporates:
     *  Sqrt: '<S72>/sqrt'
     */
    mainV03_56_B.LugV1[0] = sqrt(mainV03_56_B.u[0]) * mainV03_56_B.Output_o[2] *
      mainV03_56_B.sigma_wg;
    mainV03_56_B.LugV1[1] = sqrt(mainV03_56_B.u[1]) * mainV03_56_B.Output_o[2] *
      mainV03_56_B.MediumHighAltitudeIntensity;
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[2] == 0) {
      /* UnitDelay: '<S72>/Unit Delay' */
      mainV03_56_B.UnitDelay[0] = mainV03_56_DW.UnitDelay_DSTATE_l[0];
      mainV03_56_B.UnitDelay[1] = mainV03_56_DW.UnitDelay_DSTATE_l[1];
    }

    /* Product: '<S72>/Lug//V2' */
    mainV03_56_B.LugV2[0] = mainV03_56_B.Sum1_a[0] * mainV03_56_B.UnitDelay[0];

    /* Sum: '<S72>/Sum' */
    mainV03_56_B.Sum_da[0] = mainV03_56_B.LugV2[0] + mainV03_56_B.LugV1[0];

    /* Product: '<S72>/Lug//V2' */
    mainV03_56_B.LugV2[1] = mainV03_56_B.Sum1_a[1] * mainV03_56_B.UnitDelay[1];

    /* Sum: '<S72>/Sum' */
    mainV03_56_B.Sum_da[1] = mainV03_56_B.LugV2[1] + mainV03_56_B.LugV1[1];
  }

  /* End of Outputs for SubSystem: '<S56>/Hwgw(z)' */

  /* Outputs for Enabled SubSystem: '<S55>/Hqgw' incorporates:
   *  EnablePort: '<S68>/Enable'
   */
  if (mainV03_56_DW.Hqgw_MODE) {
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* Gain: '<S68>/dt1' */
      mainV03_56_B.dt1_k = mainV03_56_P.dt1_Gain * mainV03_56_B.UnitConversion_l;
    }

    /* Product: '<S68>/w1' */
    mainV03_56_B.aq = rtb_sincos_o2_i_idx_1 / mainV03_56_B.dt1_k;

    /* Gain: '<S68>/dt' */
    mainV03_56_B.dt_f = mainV03_56_P.dt_Gain_j * mainV03_56_B.aq;

    /* Sum: '<S68>/Sum2' incorporates:
     *  Constant: '<S68>/Constant'
     */
    mainV03_56_B.Sum2_l = mainV03_56_P.Constant_Value_k - mainV03_56_B.dt_f;
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[2] == 0) {
      /* UnitDelay: '<S68>/Unit Delay' */
      mainV03_56_B.UnitDelay_n[0] = mainV03_56_DW.UnitDelay_DSTATE_fj[0];
      mainV03_56_B.UnitDelay_n[1] = mainV03_56_DW.UnitDelay_DSTATE_fj[1];

      /* UnitDelay: '<S68>/Unit Delay1' */
      mainV03_56_B.UnitDelay1_j[0] = mainV03_56_DW.UnitDelay1_DSTATE_j[0];
      mainV03_56_B.UnitDelay1_j[1] = mainV03_56_DW.UnitDelay1_DSTATE_j[1];
    }

    /* Product: '<S68>/Lug//V2' */
    mainV03_56_B.LugV2_c3[0] = mainV03_56_B.Sum2_l * mainV03_56_B.UnitDelay_n[0];
    mainV03_56_B.LugV2_c3[1] = mainV03_56_B.Sum2_l * mainV03_56_B.UnitDelay_n[1];

    /* Sum: '<S68>/Sum3' */
    mainV03_56_B.Sum3_p[0] = mainV03_56_B.Sum_da[0] - mainV03_56_B.UnitDelay1_j
      [0];

    /* Product: '<S68>/w2' */
    mainV03_56_B.w2_b[0] = mainV03_56_B.Sum3_p[0] / mainV03_56_B.dt1_k;

    /* Sum: '<S68>/Sum1' */
    mainV03_56_B.Sum1_gx[0] = mainV03_56_B.LugV2_c3[0] + mainV03_56_B.w2_b[0];

    /* Sum: '<S68>/Sum3' */
    mainV03_56_B.Sum3_p[1] = mainV03_56_B.Sum_da[1] - mainV03_56_B.UnitDelay1_j
      [1];

    /* Product: '<S68>/w2' */
    mainV03_56_B.w2_b[1] = mainV03_56_B.Sum3_p[1] / mainV03_56_B.dt1_k;

    /* Sum: '<S68>/Sum1' */
    mainV03_56_B.Sum1_gx[1] = mainV03_56_B.LugV2_c3[1] + mainV03_56_B.w2_b[1];
  }

  /* End of Outputs for SubSystem: '<S55>/Hqgw' */

  /* Saturate: '<S74>/Limit Height h<1000ft' */
  if (rtb_sincos_o2_i_idx_0 > mainV03_56_P.LimitHeighth1000ft_UpperSat) {
    mainV03_56_B.LimitHeighth1000ft = mainV03_56_P.LimitHeighth1000ft_UpperSat;
  } else if (rtb_sincos_o2_i_idx_0 < mainV03_56_P.LimitHeighth1000ft_LowerSat) {
    mainV03_56_B.LimitHeighth1000ft = mainV03_56_P.LimitHeighth1000ft_LowerSat;
  } else {
    mainV03_56_B.LimitHeighth1000ft = rtb_sincos_o2_i_idx_0;
  }

  /* End of Saturate: '<S74>/Limit Height h<1000ft' */

  /* Fcn: '<S74>/Low Altitude Intensity' */
  rtb_IC4_idx_2 = 0.000823 * mainV03_56_B.LimitHeighth1000ft + 0.177;
  if (rtb_IC4_idx_2 < 0.0) {
    rtb_IC4_idx_2 = -rt_powd_snf(-rtb_IC4_idx_2, 0.4);
  } else {
    rtb_IC4_idx_2 = rt_powd_snf(rtb_IC4_idx_2, 0.4);
  }

  /* Product: '<S74>/sigma_ug, sigma_vg' incorporates:
   *  Fcn: '<S74>/Low Altitude Intensity'
   */
  mainV03_56_B.sigma_ugsigma_vg = 1.0 / rtb_IC4_idx_2 * mainV03_56_B.sigma_wg;

  /* Fcn: '<S91>/Low Altitude Scale Length' */
  rtb_IC4_idx_2 = 0.000823 * mainV03_56_B.LimitFunction10ftto1000ft + 0.177;
  if (rtb_IC4_idx_2 < 0.0) {
    rtb_IC4_idx_2 = -rt_powd_snf(-rtb_IC4_idx_2, 1.2);
  } else {
    rtb_IC4_idx_2 = rt_powd_snf(rtb_IC4_idx_2, 1.2);
  }

  rtb_IC4_idx_2 = mainV03_56_B.LimitFunction10ftto1000ft / rtb_IC4_idx_2;

  /* End of Fcn: '<S91>/Low Altitude Scale Length' */

  /* Gain: '<S62>/Lv' */
  mainV03_56_B.Lv[0] = mainV03_56_P.Lv_Gain * rtb_IC4_idx_2;
  mainV03_56_B.Lv[1] = mainV03_56_P.Lv_Gain * mainV03_56_B.UnitConversion_j;

  /* Outputs for Enabled SubSystem: '<S55>/Hrgw' incorporates:
   *  EnablePort: '<S69>/Enable'
   */
  /* Outputs for Enabled SubSystem: '<S56>/Hvgw(z)' incorporates:
   *  EnablePort: '<S71>/Enable'
   */
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    if (rtmIsMajorTimeStep(mainV03_56_M)) {
      /* Constant: '<S56>/Constant3' */
      if (mainV03_56_P.DrydenWindTurbulenceModelDiscreteqr_T_on > 0.0) {
        if (!mainV03_56_DW.Hvgwz_MODE) {
          /* InitializeConditions for UnitDelay: '<S71>/Unit Delay' */
          mainV03_56_DW.UnitDelay_DSTATE_o[0] =
            mainV03_56_P.UnitDelay_InitialCondition_e;
          mainV03_56_DW.UnitDelay_DSTATE_o[1] =
            mainV03_56_P.UnitDelay_InitialCondition_e;
          mainV03_56_DW.Hvgwz_MODE = true;
        }
      } else {
        if (mainV03_56_DW.Hvgwz_MODE) {
          /* Disable for Outport: '<S71>/vgw' */
          mainV03_56_B.Sum_nt[0] = mainV03_56_P.vgw_Y0;
          mainV03_56_B.Sum_nt[1] = mainV03_56_P.vgw_Y0;
          mainV03_56_DW.Hvgwz_MODE = false;
        }
      }
    }

    if (rtmIsMajorTimeStep(mainV03_56_M)) {
      /* Constant: '<S55>/Constant3' */
      if (mainV03_56_P.DrydenWindTurbulenceModelDiscreteqr_T_on > 0.0) {
        if (!mainV03_56_DW.Hrgw_MODE) {
          /* InitializeConditions for UnitDelay: '<S69>/Unit Delay' */
          mainV03_56_DW.UnitDelay_DSTATE_ju[0] =
            mainV03_56_P.UnitDelay_InitialCondition_k;

          /* InitializeConditions for UnitDelay: '<S69>/Unit Delay1' */
          mainV03_56_DW.UnitDelay1_DSTATE_ml[0] =
            mainV03_56_P.UnitDelay1_InitialCondition_h;

          /* InitializeConditions for UnitDelay: '<S69>/Unit Delay' */
          mainV03_56_DW.UnitDelay_DSTATE_ju[1] =
            mainV03_56_P.UnitDelay_InitialCondition_k;

          /* InitializeConditions for UnitDelay: '<S69>/Unit Delay1' */
          mainV03_56_DW.UnitDelay1_DSTATE_ml[1] =
            mainV03_56_P.UnitDelay1_InitialCondition_h;
          mainV03_56_DW.Hrgw_MODE = true;
        }
      } else {
        if (mainV03_56_DW.Hrgw_MODE) {
          /* Disable for Outport: '<S69>/rgw' */
          mainV03_56_B.Sum1_p2[0] = mainV03_56_P.rgw_Y0;
          mainV03_56_B.Sum1_p2[1] = mainV03_56_P.rgw_Y0;
          mainV03_56_DW.Hrgw_MODE = false;
        }
      }

      /* End of Constant: '<S55>/Constant3' */
    }
  }

  /* End of Outputs for SubSystem: '<S55>/Hrgw' */
  if (mainV03_56_DW.Hvgwz_MODE) {
    /* Product: '<S71>/V//Lvg' */
    mainV03_56_B.VLvg[0] = rtb_sincos_o2_i_idx_1 / mainV03_56_B.Lv[0];

    /* Gain: '<S71>/2' */
    mainV03_56_B.u_h[0] = mainV03_56_P.u_Gain_g * mainV03_56_B.VLvg[0];

    /* Gain: '<S71>/dt' */
    mainV03_56_B.dt_g[0] = mainV03_56_P.dt_Gain_f * mainV03_56_B.VLvg[0];

    /* Sum: '<S71>/Sum1' incorporates:
     *  Constant: '<S71>/Constant'
     */
    mainV03_56_B.Sum1_bg[0] = mainV03_56_P.Constant_Value_p - mainV03_56_B.dt_g
      [0];

    /* Product: '<S71>/V//Lvg' */
    mainV03_56_B.VLvg[1] = rtb_sincos_o2_i_idx_1 / mainV03_56_B.Lv[1];

    /* Gain: '<S71>/2' */
    mainV03_56_B.u_h[1] = mainV03_56_P.u_Gain_g * mainV03_56_B.VLvg[1];

    /* Gain: '<S71>/dt' */
    mainV03_56_B.dt_g[1] = mainV03_56_P.dt_Gain_f * mainV03_56_B.VLvg[1];

    /* Sum: '<S71>/Sum1' incorporates:
     *  Constant: '<S71>/Constant'
     */
    mainV03_56_B.Sum1_bg[1] = mainV03_56_P.Constant_Value_p - mainV03_56_B.dt_g
      [1];

    /* Product: '<S71>/Lug//V1' incorporates:
     *  Sqrt: '<S71>/sqrt'
     */
    mainV03_56_B.LugV1_m[0] = sqrt(mainV03_56_B.u_h[0]) * mainV03_56_B.Output_o
      [1] * mainV03_56_B.sigma_ugsigma_vg;
    mainV03_56_B.LugV1_m[1] = sqrt(mainV03_56_B.u_h[1]) * mainV03_56_B.Output_o
      [1] * mainV03_56_B.MediumHighAltitudeIntensity;
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[2] == 0) {
      /* UnitDelay: '<S71>/Unit Delay' */
      mainV03_56_B.UnitDelay_d[0] = mainV03_56_DW.UnitDelay_DSTATE_o[0];
      mainV03_56_B.UnitDelay_d[1] = mainV03_56_DW.UnitDelay_DSTATE_o[1];
    }

    /* Product: '<S71>/Lug//V2' */
    mainV03_56_B.LugV2_d[0] = mainV03_56_B.Sum1_bg[0] *
      mainV03_56_B.UnitDelay_d[0];

    /* Sum: '<S71>/Sum' */
    mainV03_56_B.Sum_nt[0] = mainV03_56_B.LugV2_d[0] + mainV03_56_B.LugV1_m[0];

    /* Product: '<S71>/Lug//V2' */
    mainV03_56_B.LugV2_d[1] = mainV03_56_B.Sum1_bg[1] *
      mainV03_56_B.UnitDelay_d[1];

    /* Sum: '<S71>/Sum' */
    mainV03_56_B.Sum_nt[1] = mainV03_56_B.LugV2_d[1] + mainV03_56_B.LugV1_m[1];
  }

  /* End of Outputs for SubSystem: '<S56>/Hvgw(z)' */

  /* Outputs for Enabled SubSystem: '<S55>/Hrgw' incorporates:
   *  EnablePort: '<S69>/Enable'
   */
  if (mainV03_56_DW.Hrgw_MODE) {
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* Gain: '<S69>/dt1' */
      mainV03_56_B.dt1 = mainV03_56_P.dt1_Gain_l * mainV03_56_B.UnitConversion_l;
    }

    /* Product: '<S69>/w1' */
    mainV03_56_B.ar = rtb_sincos_o2_i_idx_1 / mainV03_56_B.dt1;

    /* Gain: '<S69>/dt' */
    mainV03_56_B.dt_i = mainV03_56_P.dt_Gain_k * mainV03_56_B.ar;

    /* Sum: '<S69>/Sum2' incorporates:
     *  Constant: '<S69>/Constant'
     */
    mainV03_56_B.Sum2_fh = mainV03_56_P.Constant_Value_l - mainV03_56_B.dt_i;
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[2] == 0) {
      /* UnitDelay: '<S69>/Unit Delay' */
      mainV03_56_B.UnitDelay_b[0] = mainV03_56_DW.UnitDelay_DSTATE_ju[0];
      mainV03_56_B.UnitDelay_b[1] = mainV03_56_DW.UnitDelay_DSTATE_ju[1];

      /* UnitDelay: '<S69>/Unit Delay1' */
      mainV03_56_B.UnitDelay1[0] = mainV03_56_DW.UnitDelay1_DSTATE_ml[0];
      mainV03_56_B.UnitDelay1[1] = mainV03_56_DW.UnitDelay1_DSTATE_ml[1];
    }

    /* Product: '<S69>/Lug//V2' */
    mainV03_56_B.LugV2_p[0] = mainV03_56_B.Sum2_fh * mainV03_56_B.UnitDelay_b[0];
    mainV03_56_B.LugV2_p[1] = mainV03_56_B.Sum2_fh * mainV03_56_B.UnitDelay_b[1];

    /* Sum: '<S69>/Sum3' */
    mainV03_56_B.Sum3_jj[0] = mainV03_56_B.Sum_nt[0] - mainV03_56_B.UnitDelay1[0];

    /* Product: '<S69>/w2' */
    mainV03_56_B.w2[0] = mainV03_56_B.Sum3_jj[0] / mainV03_56_B.dt1;

    /* Sum: '<S69>/Sum1' */
    mainV03_56_B.Sum1_p2[0] = mainV03_56_B.LugV2_p[0] + mainV03_56_B.w2[0];

    /* Sum: '<S69>/Sum3' */
    mainV03_56_B.Sum3_jj[1] = mainV03_56_B.Sum_nt[1] - mainV03_56_B.UnitDelay1[1];

    /* Product: '<S69>/w2' */
    mainV03_56_B.w2[1] = mainV03_56_B.Sum3_jj[1] / mainV03_56_B.dt1;

    /* Sum: '<S69>/Sum1' */
    mainV03_56_B.Sum1_p2[1] = mainV03_56_B.LugV2_p[1] + mainV03_56_B.w2[1];
  }

  /* End of Outputs for SubSystem: '<S55>/Hrgw' */

  /* Outputs for Enabled SubSystem: '<S56>/Hugw(z)' incorporates:
   *  EnablePort: '<S70>/Enable'
   */
  if ((rtmIsMajorTimeStep(mainV03_56_M) &&
       mainV03_56_M->Timing.TaskCounters.TID[1] == 0) && rtmIsMajorTimeStep
      (mainV03_56_M)) {
    /* Constant: '<S56>/Constant3' */
    if (mainV03_56_P.DrydenWindTurbulenceModelDiscreteqr_T_on > 0.0) {
      if (!mainV03_56_DW.Hugwz_MODE) {
        /* InitializeConditions for UnitDelay: '<S70>/Unit Delay' */
        mainV03_56_DW.UnitDelay_DSTATE_p[0] =
          mainV03_56_P.UnitDelay_InitialCondition_m;
        mainV03_56_DW.UnitDelay_DSTATE_p[1] =
          mainV03_56_P.UnitDelay_InitialCondition_m;
        mainV03_56_DW.Hugwz_MODE = true;
      }
    } else {
      if (mainV03_56_DW.Hugwz_MODE) {
        /* Disable for Outport: '<S70>/ugw' */
        mainV03_56_B.Sum_eq[0] = mainV03_56_P.ugw_Y0;
        mainV03_56_B.Sum_eq[1] = mainV03_56_P.ugw_Y0;
        mainV03_56_DW.Hugwz_MODE = false;
      }
    }
  }

  if (mainV03_56_DW.Hugwz_MODE) {
    /* Product: '<S70>/V//Lug' */
    mainV03_56_B.VLug[0] = rtb_sincos_o2_i_idx_1 / rtb_IC4_idx_2;
    mainV03_56_B.VLug[1] = rtb_sincos_o2_i_idx_1 / mainV03_56_B.UnitConversion_j;

    /* Gain: '<S70>/2' */
    mainV03_56_B.u_i[0] = mainV03_56_P.u_Gain_a * mainV03_56_B.VLug[0];

    /* Gain: '<S70>/dt' */
    mainV03_56_B.dt_l[0] = mainV03_56_P.dt_Gain_jf * mainV03_56_B.VLug[0];

    /* Sum: '<S70>/Sum1' incorporates:
     *  Constant: '<S70>/Constant'
     */
    mainV03_56_B.Sum1_o5[0] = mainV03_56_P.Constant_Value_b - mainV03_56_B.dt_l
      [0];

    /* Gain: '<S70>/2' */
    mainV03_56_B.u_i[1] = mainV03_56_P.u_Gain_a * mainV03_56_B.VLug[1];

    /* Gain: '<S70>/dt' */
    mainV03_56_B.dt_l[1] = mainV03_56_P.dt_Gain_jf * mainV03_56_B.VLug[1];

    /* Sum: '<S70>/Sum1' incorporates:
     *  Constant: '<S70>/Constant'
     */
    mainV03_56_B.Sum1_o5[1] = mainV03_56_P.Constant_Value_b - mainV03_56_B.dt_l
      [1];

    /* Product: '<S70>/Lug//V1' incorporates:
     *  Sqrt: '<S70>/sqrt'
     */
    mainV03_56_B.LugV1_l[0] = sqrt(mainV03_56_B.u_i[0]) * mainV03_56_B.Output_o
      [0] * mainV03_56_B.sigma_ugsigma_vg;
    mainV03_56_B.LugV1_l[1] = sqrt(mainV03_56_B.u_i[1]) * mainV03_56_B.Output_o
      [0] * mainV03_56_B.MediumHighAltitudeIntensity;
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[2] == 0) {
      /* UnitDelay: '<S70>/Unit Delay' */
      mainV03_56_B.UnitDelay_k[0] = mainV03_56_DW.UnitDelay_DSTATE_p[0];
      mainV03_56_B.UnitDelay_k[1] = mainV03_56_DW.UnitDelay_DSTATE_p[1];
    }

    /* Product: '<S70>/Lug//V2' */
    mainV03_56_B.LugV2_c[0] = mainV03_56_B.Sum1_o5[0] *
      mainV03_56_B.UnitDelay_k[0];

    /* Sum: '<S70>/Sum' */
    mainV03_56_B.Sum_eq[0] = mainV03_56_B.LugV2_c[0] + mainV03_56_B.LugV1_l[0];

    /* Product: '<S70>/Lug//V2' */
    mainV03_56_B.LugV2_c[1] = mainV03_56_B.Sum1_o5[1] *
      mainV03_56_B.UnitDelay_k[1];

    /* Sum: '<S70>/Sum' */
    mainV03_56_B.Sum_eq[1] = mainV03_56_B.LugV2_c[1] + mainV03_56_B.LugV1_l[1];
  }

  /* End of Outputs for SubSystem: '<S56>/Hugw(z)' */

  /* If: '<S60>/if Height < Max low altitude  elseif Height > Min isotropic altitude ' */
  if (rtmIsMajorTimeStep(mainV03_56_M)) {
    if (rtb_sincos_o2_i_idx_0 <= 1000.0) {
      rtAction = 0;
    } else if (rtb_sincos_o2_i_idx_0 >= 2000.0) {
      rtAction = 1;
    } else {
      rtAction = 2;
    }

    mainV03_56_DW.ifHeightMaxlowaltitudeelseifHeightMinisotropicaltitude_ActiveSubsystem
      = rtAction;
  } else {
    rtAction =
      mainV03_56_DW.ifHeightMaxlowaltitudeelseifHeightMinisotropicaltitude_ActiveSubsystem;
  }

  switch (rtAction) {
   case 0:
    /* Outputs for IfAction SubSystem: '<S60>/Low altitude  rates' incorporates:
     *  ActionPort: '<S76>/Action Port'
     */
    /* SignalConversion: '<S81>/ConcatBufferAtVector ConcatenateIn3' */
    mainV03_56_B.VectorConcatenate_aw[2] = mainV03_56_B.Sum1_p2[0];

    /* Trigonometry: '<S82>/Trigonometric Function1' */
    rtb_sincos_o2_i_idx_1 = sin(mainV03_56_B.UnitConversion);
    rtb_IC4_idx_2 = cos(mainV03_56_B.UnitConversion);

    /* Product: '<S82>/Product2' */
    mainV03_56_B.Product2_id[0] = mainV03_56_B.Sum_ax[0] * rtb_IC4_idx_2;
    mainV03_56_B.Product2_id[1] = mainV03_56_B.Sum1_gx[0] * rtb_IC4_idx_2;

    /* Product: '<S82>/Product1' */
    mainV03_56_B.Product1_dl[0] = rtb_sincos_o2_i_idx_1 * mainV03_56_B.Sum_ax[0];
    mainV03_56_B.Product1_dl[1] = rtb_sincos_o2_i_idx_1 * mainV03_56_B.Sum1_gx[0];

    /* Sum: '<S82>/Sum' */
    mainV03_56_B.VectorConcatenate_aw[0] = mainV03_56_B.Product2_id[0] -
      mainV03_56_B.Product1_dl[1];

    /* Sum: '<S82>/Sum1' */
    mainV03_56_B.VectorConcatenate_aw[1] = mainV03_56_B.Product2_id[1] +
      mainV03_56_B.Product1_dl[0];
    for (s155_iter = 0; s155_iter < 3; s155_iter++) {
      /* Product: '<S81>/Product' */
      mainV03_56_B.Product_kg[s155_iter] = 0.0;
      mainV03_56_B.Product_kg[s155_iter] +=
        mainV03_56_B.plantData.DCM_body_earth[s155_iter] *
        mainV03_56_B.VectorConcatenate_aw[0];
      mainV03_56_B.Product_kg[s155_iter] +=
        mainV03_56_B.plantData.DCM_body_earth[s155_iter + 3] *
        mainV03_56_B.VectorConcatenate_aw[1];
      mainV03_56_B.Product_kg[s155_iter] +=
        mainV03_56_B.plantData.DCM_body_earth[s155_iter + 6] *
        mainV03_56_B.VectorConcatenate_aw[2];

      /* Reshape: '<S81>/Reshape1' */
      mainV03_56_B.Merge_n[s155_iter] = mainV03_56_B.Product_kg[s155_iter];
    }

    /* End of Outputs for SubSystem: '<S60>/Low altitude  rates' */
    break;

   case 1:
    /* Outputs for IfAction SubSystem: '<S60>/Medium//High  altitude rates' incorporates:
     *  ActionPort: '<S77>/Action Port'
     */
    /* Gain: '<S77>/Gain' */
    mainV03_56_B.Merge_n[0] = mainV03_56_P.Gain_Gain * mainV03_56_B.Sum_ax[1];
    mainV03_56_B.Merge_n[1] = mainV03_56_P.Gain_Gain * mainV03_56_B.Sum1_gx[1];
    mainV03_56_B.Merge_n[2] = mainV03_56_P.Gain_Gain * mainV03_56_B.Sum1_p2[1];

    /* End of Outputs for SubSystem: '<S60>/Medium//High  altitude rates' */
    break;

   case 2:
    /* Outputs for IfAction SubSystem: '<S60>/Interpolate  rates' incorporates:
     *  ActionPort: '<S75>/Action Port'
     */
    /* Trigonometry: '<S80>/Trigonometric Function' */
    rtb_sincos_o2_i_idx_1 = sin(mainV03_56_B.UnitConversion);
    rtb_IC4_idx_2 = cos(mainV03_56_B.UnitConversion);

    /* Product: '<S80>/Product2' */
    mainV03_56_B.Product2_my[0] = mainV03_56_B.Sum_ax[0] * rtb_IC4_idx_2;
    mainV03_56_B.Product2_my[1] = mainV03_56_B.Sum1_gx[0] * rtb_IC4_idx_2;

    /* Product: '<S80>/Product1' */
    mainV03_56_B.Product1_k0[0] = rtb_sincos_o2_i_idx_1 * mainV03_56_B.Sum_ax[0];
    mainV03_56_B.Product1_k0[1] = rtb_sincos_o2_i_idx_1 * mainV03_56_B.Sum1_gx[0];

    /* Sum: '<S80>/Sum' */
    mainV03_56_B.VectorConcatenate_p[0] = mainV03_56_B.Product2_my[0] -
      mainV03_56_B.Product1_k0[1];

    /* Sum: '<S80>/Sum1' */
    mainV03_56_B.VectorConcatenate_p[1] = mainV03_56_B.Product2_my[1] +
      mainV03_56_B.Product1_k0[0];

    /* SignalConversion: '<S79>/ConcatBufferAtVector ConcatenateIn3' */
    mainV03_56_B.VectorConcatenate_p[2] = mainV03_56_B.Sum1_p2[0];

    /* Product: '<S79>/Product' */
    for (rtb_Sum1_eh = 0; rtb_Sum1_eh < 3; rtb_Sum1_eh++) {
      mainV03_56_B.Product_gv[rtb_Sum1_eh] = 0.0;
      mainV03_56_B.Product_gv[rtb_Sum1_eh] +=
        mainV03_56_B.plantData.DCM_body_earth[rtb_Sum1_eh] *
        mainV03_56_B.VectorConcatenate_p[0];
      mainV03_56_B.Product_gv[rtb_Sum1_eh] +=
        mainV03_56_B.plantData.DCM_body_earth[rtb_Sum1_eh + 3] *
        mainV03_56_B.VectorConcatenate_p[1];
      mainV03_56_B.Product_gv[rtb_Sum1_eh] +=
        mainV03_56_B.plantData.DCM_body_earth[rtb_Sum1_eh + 6] *
        mainV03_56_B.VectorConcatenate_p[2];
    }

    /* End of Product: '<S79>/Product' */

    /* Sum: '<S75>/Sum2' */
    mainV03_56_B.Sum2_g5[0] = mainV03_56_B.Sum_ax[1] - mainV03_56_B.Product_gv[0];
    mainV03_56_B.Sum2_g5[1] = mainV03_56_B.Sum1_gx[1] - mainV03_56_B.Product_gv
      [1];
    mainV03_56_B.Sum2_g5[2] = mainV03_56_B.Sum1_p2[1] - mainV03_56_B.Product_gv
      [2];

    /* Sum: '<S75>/Sum1' incorporates:
     *  Constant: '<S75>/max_height_low'
     */
    mainV03_56_B.Sum1_jk = rtb_sincos_o2_i_idx_0 -
      mainV03_56_P.max_height_low_Value;
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* Sum: '<S75>/Sum' incorporates:
       *  Constant: '<S75>/max_height_low'
       *  Constant: '<S75>/min_height_high'
       */
      mainV03_56_B.Sum_eg = mainV03_56_P.min_height_high_Value -
        mainV03_56_P.max_height_low_Value;
    }

    /* Product: '<S75>/Product1' */
    mainV03_56_B.Product1_h5[0] = mainV03_56_B.Sum2_g5[0] * mainV03_56_B.Sum1_jk
      / mainV03_56_B.Sum_eg;

    /* Sum: '<S75>/Sum3' */
    mainV03_56_B.Merge_n[0] = mainV03_56_B.Product_gv[0] +
      mainV03_56_B.Product1_h5[0];

    /* Product: '<S75>/Product1' */
    mainV03_56_B.Product1_h5[1] = mainV03_56_B.Sum2_g5[1] * mainV03_56_B.Sum1_jk
      / mainV03_56_B.Sum_eg;

    /* Sum: '<S75>/Sum3' */
    mainV03_56_B.Merge_n[1] = mainV03_56_B.Product_gv[1] +
      mainV03_56_B.Product1_h5[1];

    /* Product: '<S75>/Product1' */
    mainV03_56_B.Product1_h5[2] = mainV03_56_B.Sum2_g5[2] * mainV03_56_B.Sum1_jk
      / mainV03_56_B.Sum_eg;

    /* Sum: '<S75>/Sum3' */
    mainV03_56_B.Merge_n[2] = mainV03_56_B.Product_gv[2] +
      mainV03_56_B.Product1_h5[2];

    /* End of Outputs for SubSystem: '<S60>/Interpolate  rates' */
    break;
  }

  /* End of If: '<S60>/if Height < Max low altitude  elseif Height > Min isotropic altitude ' */

  /* If: '<S61>/if Height < Max low altitude  elseif Height > Min isotropic altitude ' */
  if (rtmIsMajorTimeStep(mainV03_56_M)) {
    if (rtb_sincos_o2_i_idx_0 <= 1000.0) {
      rtAction = 0;
    } else if (rtb_sincos_o2_i_idx_0 >= 2000.0) {
      rtAction = 1;
    } else {
      rtAction = 2;
    }

    mainV03_56_DW.ifHeightMaxlowaltitudeelseifHeightMinisotropicaltitude_ActiveSubsystem_d
      = rtAction;
  } else {
    rtAction =
      mainV03_56_DW.ifHeightMaxlowaltitudeelseifHeightMinisotropicaltitude_ActiveSubsystem_d;
  }

  switch (rtAction) {
   case 0:
    /* Outputs for IfAction SubSystem: '<S61>/Low altitude  velocities' incorporates:
     *  ActionPort: '<S84>/Action Port'
     */
    /* SignalConversion: '<S89>/ConcatBufferAtVector ConcatenateIn3' */
    mainV03_56_B.VectorConcatenate_a[2] = mainV03_56_B.Sum_da[0];

    /* Trigonometry: '<S90>/Trigonometric Function' */
    rtb_sincos_o2_i_idx_0 = sin(mainV03_56_B.UnitConversion);
    rtb_sincos_o2_i_idx_1 = cos(mainV03_56_B.UnitConversion);

    /* Product: '<S90>/Product2' */
    mainV03_56_B.Product2_oo[0] = mainV03_56_B.Sum_eq[0] * rtb_sincos_o2_i_idx_1;
    mainV03_56_B.Product2_oo[1] = mainV03_56_B.Sum_nt[0] * rtb_sincos_o2_i_idx_1;

    /* Product: '<S90>/Product1' */
    mainV03_56_B.Product1_o0[0] = rtb_sincos_o2_i_idx_0 * mainV03_56_B.Sum_eq[0];
    mainV03_56_B.Product1_o0[1] = rtb_sincos_o2_i_idx_0 * mainV03_56_B.Sum_nt[0];

    /* Sum: '<S90>/Sum' */
    mainV03_56_B.VectorConcatenate_a[0] = mainV03_56_B.Product2_oo[0] -
      mainV03_56_B.Product1_o0[1];

    /* Sum: '<S90>/Sum1' */
    mainV03_56_B.VectorConcatenate_a[1] = mainV03_56_B.Product2_oo[1] +
      mainV03_56_B.Product1_o0[0];
    for (s155_iter = 0; s155_iter < 3; s155_iter++) {
      /* Product: '<S89>/Product' */
      mainV03_56_B.Product_mb[s155_iter] = 0.0;
      mainV03_56_B.Product_mb[s155_iter] +=
        mainV03_56_B.plantData.DCM_body_earth[s155_iter] *
        mainV03_56_B.VectorConcatenate_a[0];
      mainV03_56_B.Product_mb[s155_iter] +=
        mainV03_56_B.plantData.DCM_body_earth[s155_iter + 3] *
        mainV03_56_B.VectorConcatenate_a[1];
      mainV03_56_B.Product_mb[s155_iter] +=
        mainV03_56_B.plantData.DCM_body_earth[s155_iter + 6] *
        mainV03_56_B.VectorConcatenate_a[2];

      /* Reshape: '<S89>/Reshape1' */
      mainV03_56_B.Merge_i[s155_iter] = mainV03_56_B.Product_mb[s155_iter];
    }

    /* End of Outputs for SubSystem: '<S61>/Low altitude  velocities' */
    break;

   case 1:
    /* Outputs for IfAction SubSystem: '<S61>/Medium//High  altitude velocities' incorporates:
     *  ActionPort: '<S85>/Action Port'
     */
    /* Gain: '<S85>/Gain' */
    mainV03_56_B.Merge_i[0] = mainV03_56_P.Gain_Gain_c * mainV03_56_B.Sum_eq[1];
    mainV03_56_B.Merge_i[1] = mainV03_56_P.Gain_Gain_c * mainV03_56_B.Sum_nt[1];
    mainV03_56_B.Merge_i[2] = mainV03_56_P.Gain_Gain_c * mainV03_56_B.Sum_da[1];

    /* End of Outputs for SubSystem: '<S61>/Medium//High  altitude velocities' */
    break;

   case 2:
    /* Outputs for IfAction SubSystem: '<S61>/Interpolate  velocities' incorporates:
     *  ActionPort: '<S83>/Action Port'
     */
    /* Trigonometry: '<S88>/Trigonometric Function' */
    rtb_sincos_o2_i_idx_1 = sin(mainV03_56_B.UnitConversion);
    rtb_IC4_idx_2 = cos(mainV03_56_B.UnitConversion);

    /* Product: '<S88>/Product2' */
    mainV03_56_B.Product2_i[0] = mainV03_56_B.Sum_eq[0] * rtb_IC4_idx_2;
    mainV03_56_B.Product2_i[1] = mainV03_56_B.Sum_nt[0] * rtb_IC4_idx_2;

    /* Product: '<S88>/Product1' */
    mainV03_56_B.Product1_mk[0] = rtb_sincos_o2_i_idx_1 * mainV03_56_B.Sum_eq[0];
    mainV03_56_B.Product1_mk[1] = rtb_sincos_o2_i_idx_1 * mainV03_56_B.Sum_nt[0];

    /* Sum: '<S88>/Sum' */
    mainV03_56_B.VectorConcatenate_j3[0] = mainV03_56_B.Product2_i[0] -
      mainV03_56_B.Product1_mk[1];

    /* Sum: '<S88>/Sum1' */
    mainV03_56_B.VectorConcatenate_j3[1] = mainV03_56_B.Product2_i[1] +
      mainV03_56_B.Product1_mk[0];

    /* SignalConversion: '<S87>/ConcatBufferAtVector ConcatenateIn3' */
    mainV03_56_B.VectorConcatenate_j3[2] = mainV03_56_B.Sum_da[0];

    /* Product: '<S87>/Product' */
    for (rtb_Sum1_eh = 0; rtb_Sum1_eh < 3; rtb_Sum1_eh++) {
      mainV03_56_B.Product_lt[rtb_Sum1_eh] = 0.0;
      mainV03_56_B.Product_lt[rtb_Sum1_eh] +=
        mainV03_56_B.plantData.DCM_body_earth[rtb_Sum1_eh] *
        mainV03_56_B.VectorConcatenate_j3[0];
      mainV03_56_B.Product_lt[rtb_Sum1_eh] +=
        mainV03_56_B.plantData.DCM_body_earth[rtb_Sum1_eh + 3] *
        mainV03_56_B.VectorConcatenate_j3[1];
      mainV03_56_B.Product_lt[rtb_Sum1_eh] +=
        mainV03_56_B.plantData.DCM_body_earth[rtb_Sum1_eh + 6] *
        mainV03_56_B.VectorConcatenate_j3[2];
    }

    /* End of Product: '<S87>/Product' */

    /* Sum: '<S83>/Sum2' */
    mainV03_56_B.Sum2_g[0] = mainV03_56_B.Sum_eq[1] - mainV03_56_B.Product_lt[0];
    mainV03_56_B.Sum2_g[1] = mainV03_56_B.Sum_nt[1] - mainV03_56_B.Product_lt[1];
    mainV03_56_B.Sum2_g[2] = mainV03_56_B.Sum_da[1] - mainV03_56_B.Product_lt[2];

    /* Sum: '<S83>/Sum1' incorporates:
     *  Constant: '<S83>/max_height_low'
     */
    mainV03_56_B.Sum1_jc = rtb_sincos_o2_i_idx_0 -
      mainV03_56_P.max_height_low_Value_i;
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* Sum: '<S83>/Sum' incorporates:
       *  Constant: '<S83>/max_height_low'
       *  Constant: '<S83>/min_height_high'
       */
      mainV03_56_B.Sum_em = mainV03_56_P.min_height_high_Value_n -
        mainV03_56_P.max_height_low_Value_i;
    }

    /* Product: '<S83>/Product1' */
    mainV03_56_B.Product1_hj[0] = mainV03_56_B.Sum2_g[0] * mainV03_56_B.Sum1_jc /
      mainV03_56_B.Sum_em;

    /* Sum: '<S83>/Sum3' */
    mainV03_56_B.Merge_i[0] = mainV03_56_B.Product_lt[0] +
      mainV03_56_B.Product1_hj[0];

    /* Product: '<S83>/Product1' */
    mainV03_56_B.Product1_hj[1] = mainV03_56_B.Sum2_g[1] * mainV03_56_B.Sum1_jc /
      mainV03_56_B.Sum_em;

    /* Sum: '<S83>/Sum3' */
    mainV03_56_B.Merge_i[1] = mainV03_56_B.Product_lt[1] +
      mainV03_56_B.Product1_hj[1];

    /* Product: '<S83>/Product1' */
    mainV03_56_B.Product1_hj[2] = mainV03_56_B.Sum2_g[2] * mainV03_56_B.Sum1_jc /
      mainV03_56_B.Sum_em;

    /* Sum: '<S83>/Sum3' */
    mainV03_56_B.Merge_i[2] = mainV03_56_B.Product_lt[2] +
      mainV03_56_B.Product1_hj[2];

    /* End of Outputs for SubSystem: '<S61>/Interpolate  velocities' */
    break;
  }

  /* End of If: '<S61>/if Height < Max low altitude  elseif Height > Min isotropic altitude ' */
  /* Unit Conversion - from: ft/s to: m/s
     Expression: output = (0.3048*input) + (0) */
  /* Unit Conversion - from: m to: ft
     Expression: output = (3.28084*input) + (0) */
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Math: '<S49>/ln(ref_height//z0)' incorporates:
     *  Constant: '<S49>/ref_height//z0'
     *
     * About '<S49>/ln(ref_height//z0)':
     *  Operator: log
     */
    mainV03_56_B.lnref_heightz0 = log(mainV03_56_P.ref_heightz0_Value);

    /* UnitConversion: '<S94>/Unit Conversion' incorporates:
     *  Constant: '<S49>/Wind Direction'
     */
    /* Unit Conversion - from: deg to: rad
       Expression: output = (0.0174533*input) + (0) */
    rtb_Switch = 0.017453292519943295 * mainV03_56_P.WindShearModel_Wdeg;

    /* Trigonometry: '<S49>/SinCos' */
    rtb_Abs1 = cos(rtb_Switch);
    rtb_Switch = sin(rtb_Switch);

    /* Gain: '<S49>/Wind speed at reference height' incorporates:
     *  Constant: '<S49>/Wdeg1'
     *  Trigonometry: '<S49>/SinCos'
     */
    mainV03_56_B.Windspeedatreferenceheight[0] =
      -mainV03_56_P.WindShearModel_W_20 * rtb_Abs1;
    mainV03_56_B.Windspeedatreferenceheight[1] =
      -mainV03_56_P.WindShearModel_W_20 * rtb_Switch;
    mainV03_56_B.Windspeedatreferenceheight[2] =
      -mainV03_56_P.WindShearModel_W_20 * mainV03_56_P.Wdeg1_Value;
  }

  /* ManualSwitch: '<S15>/Manual Switch' incorporates:
   *  Constant: '<S15>/Constant2'
   */
  if (mainV03_56_P.ManualSwitch_CurrentSetting_o == 1) {
    mainV03_56_B.envData.windVelocity[0] = mainV03_56_P.Constant2_Value_cd[0];
    mainV03_56_B.envData.windVelocity[1] = mainV03_56_P.Constant2_Value_cd[1];
    mainV03_56_B.envData.windVelocity[2] = mainV03_56_P.Constant2_Value_cd[2];
  } else {
    /* UnitConversion: '<S95>/Unit Conversion' */
    rtb_Sum1_no = 3.280839895013123 * mainV03_56_B.plantData.LLA.Altitude_m;

    /* Saturate: '<S49>/3ft-->inf' */
    if (rtb_Sum1_no > mainV03_56_P.uftinf_UpperSat) {
      mainV03_56_B.uftinf = mainV03_56_P.uftinf_UpperSat;
    } else if (rtb_Sum1_no < mainV03_56_P.uftinf_LowerSat) {
      mainV03_56_B.uftinf = mainV03_56_P.uftinf_LowerSat;
    } else {
      mainV03_56_B.uftinf = rtb_Sum1_no;
    }

    /* End of Saturate: '<S49>/3ft-->inf' */

    /* Gain: '<S49>/h//z0' */
    mainV03_56_B.hz0 = mainV03_56_P.hz0_Gain * mainV03_56_B.uftinf;

    /* Product: '<S49>/Product' incorporates:
     *  Math: '<S49>/ln(h//z0)'
     *
     * About '<S49>/ln(h//z0)':
     *  Operator: log
     */
    mainV03_56_B.Product_dl = log(mainV03_56_B.hz0) /
      mainV03_56_B.lnref_heightz0;

    /* Product: '<S49>/Product1' */
    mainV03_56_B.Product1_ek[0] = mainV03_56_B.Product_dl *
      mainV03_56_B.Windspeedatreferenceheight[0];
    mainV03_56_B.Product1_ek[1] = mainV03_56_B.Product_dl *
      mainV03_56_B.Windspeedatreferenceheight[1];
    mainV03_56_B.Product1_ek[2] = mainV03_56_B.Product_dl *
      mainV03_56_B.Windspeedatreferenceheight[2];
    for (s155_iter = 0; s155_iter < 3; s155_iter++) {
      /* Product: '<S49>/Transform from Inertial to Body axes' */
      mainV03_56_B.TransformfromInertialtoBodyaxes[s155_iter] = 0.0;
      mainV03_56_B.TransformfromInertialtoBodyaxes[s155_iter] +=
        mainV03_56_B.plantData.DCM_body_earth[s155_iter] *
        mainV03_56_B.Product1_ek[0];
      mainV03_56_B.TransformfromInertialtoBodyaxes[s155_iter] +=
        mainV03_56_B.plantData.DCM_body_earth[s155_iter + 3] *
        mainV03_56_B.Product1_ek[1];
      mainV03_56_B.TransformfromInertialtoBodyaxes[s155_iter] +=
        mainV03_56_B.plantData.DCM_body_earth[s155_iter + 6] *
        mainV03_56_B.Product1_ek[2];

      /* Sum: '<S15>/Sum' incorporates:
       *  UnitConversion: '<S48>/Unit Conversion'
       */
      mainV03_56_B.WindVelocity[s155_iter] = (0.3048 *
        mainV03_56_B.Merge_i[s155_iter] +
        mainV03_56_B.TransformfromInertialtoBodyaxes[s155_iter]) +
        mainV03_56_B.Gustmagnitude20[s155_iter];
      mainV03_56_B.envData.windVelocity[s155_iter] =
        mainV03_56_B.WindVelocity[s155_iter];
    }
  }

  /* End of ManualSwitch: '<S15>/Manual Switch' */

  /* ManualSwitch: '<S15>/Manual Switch1' incorporates:
   *  Constant: '<S15>/Constant'
   */
  if (mainV03_56_P.ManualSwitch1_CurrentSetting_g == 1) {
    mainV03_56_B.envData.windOmega[0] = mainV03_56_B.Merge_n[0];
    mainV03_56_B.envData.windOmega[1] = mainV03_56_B.Merge_n[1];
    mainV03_56_B.envData.windOmega[2] = mainV03_56_B.Merge_n[2];
  } else {
    mainV03_56_B.envData.windOmega[0] = mainV03_56_P.Constant_Value_bp[0];
    mainV03_56_B.envData.windOmega[1] = mainV03_56_P.Constant_Value_bp[1];
    mainV03_56_B.envData.windOmega[2] = mainV03_56_P.Constant_Value_bp[2];
  }

  /* End of ManualSwitch: '<S15>/Manual Switch1' */

  /* Logic: '<S96>/conjunction' incorporates:
   *  Constant: '<S96>/max_val'
   *  Constant: '<S96>/min_val'
   *  RelationalOperator: '<S96>/max_relop'
   *  RelationalOperator: '<S96>/min_relop'
   */
  mainV03_56_B.conjunction = ((mainV03_56_P.CheckAltitude_min <=
    mainV03_56_B.plantData.LLA.Altitude_m) &&
    (mainV03_56_B.plantData.LLA.Altitude_m <= mainV03_56_P.CheckAltitude_max));
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Assertion: '<S96>/Assertion' */
    utAssert(mainV03_56_B.conjunction);
  }

  /* Logic: '<S97>/conjunction' incorporates:
   *  Constant: '<S97>/max_val'
   *  Constant: '<S97>/min_val'
   *  RelationalOperator: '<S97>/max_relop'
   *  RelationalOperator: '<S97>/min_relop'
   */
  mainV03_56_B.conjunction_g = ((mainV03_56_P.CheckLatitude_min <=
    mainV03_56_B.plantData.LLA.Latitude_deg) &&
    (mainV03_56_B.plantData.LLA.Latitude_deg <= mainV03_56_P.CheckLatitude_max));
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Assertion: '<S97>/Assertion' */
    utAssert(mainV03_56_B.conjunction_g);
  }

  /* Logic: '<S98>/conjunction' incorporates:
   *  Constant: '<S98>/max_val'
   *  Constant: '<S98>/min_val'
   *  RelationalOperator: '<S98>/max_relop'
   *  RelationalOperator: '<S98>/min_relop'
   */
  mainV03_56_B.conjunction_d = ((mainV03_56_P.CheckLongitude_min <=
    mainV03_56_B.plantData.LLA.Longitude_deg) &&
    (mainV03_56_B.plantData.LLA.Longitude_deg <= mainV03_56_P.CheckLongitude_max));
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Assertion: '<S98>/Assertion' */
    utAssert(mainV03_56_B.conjunction_d);
  }

  /* Logic: '<S100>/conjunction' incorporates:
   *  Constant: '<S100>/max_val'
   *  Constant: '<S100>/min_val'
   *  RelationalOperator: '<S100>/max_relop'
   *  RelationalOperator: '<S100>/min_relop'
   */
  mainV03_56_B.conjunction_k = ((mainV03_56_P.Istimewithinmodellimits_min <=
    mainV03_56_B.Sum_n) && (mainV03_56_B.Sum_n <=
    mainV03_56_P.Istimewithinmodellimits_max));
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Assertion: '<S100>/Assertion' */
    utAssert(mainV03_56_B.conjunction_k);
  }

  /* BusCreator: '<S2>/envDataBus' */
  mainV03_56_B.envData.h_ground = 0.0;
  mainV03_56_B.envData.T = mainV03_56_B.SFunction_o1;
  mainV03_56_B.envData.a = mainV03_56_B.SFunction_o2;
  mainV03_56_B.envData.P = mainV03_56_B.SFunction_o3;
  mainV03_56_B.envData.AirDensity = mainV03_56_B.SFunction_o4;
  mainV03_56_B.envData.Gravity_ned[0] = mainV03_56_B.MatrixConcatenate[0];
  mainV03_56_B.envData.Gravity_ned[1] = mainV03_56_B.MatrixConcatenate[1];
  mainV03_56_B.envData.Gravity_ned[2] = mainV03_56_B.MatrixConcatenate[2];

  /* FromWorkspace: '<S463>/FromWs' */
  {
    real_T *pDataValues = (real_T *) mainV03_56_DW.FromWs_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *) mainV03_56_DW.FromWs_PWORK.TimePtr;
    int_T currTimeIndex = mainV03_56_DW.FromWs_IWORK.PrevIndex;
    real_T t = mainV03_56_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[13]) {
      currTimeIndex = 12;
    } else {
      if (t < pTimeValues[currTimeIndex]) {
        while (t < pTimeValues[currTimeIndex]) {
          currTimeIndex--;
        }
      } else {
        while (t >= pTimeValues[currTimeIndex + 1]) {
          currTimeIndex++;
        }
      }
    }

    mainV03_56_DW.FromWs_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 6; ++elIdx) {
              (&rtb_FromWs[0])[elIdx] = pDataValues[currTimeIndex];
              pDataValues += 14;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 6; ++elIdx) {
              (&rtb_FromWs[0])[elIdx] = pDataValues[currTimeIndex + 1];
              pDataValues += 14;
            }
          }
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;

        {
          int_T elIdx;
          for (elIdx = 0; elIdx < 6; ++elIdx) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            (&rtb_FromWs[0])[elIdx] = (real_T) rtInterpolate(d1, d2, f1, f2);
            pDataValues += 14;
          }
        }
      }
    }
  }

  /* BusCreator: '<S5>/Bus Creator' incorporates:
   *  Constant: '<S5>/Comanded position in Earth axis'
   */
  mainV03_56_B.BusCreator.roll_cmd = rtb_FromWs[0];
  mainV03_56_B.BusCreator.pitch_cmd = rtb_FromWs[1];
  mainV03_56_B.BusCreator.yaw_cmd = rtb_FromWs[2];
  mainV03_56_B.BusCreator.altitude_cmd = rtb_FromWs[3];
  mainV03_56_B.BusCreator.Throttle1_cmd = rtb_FromWs[4];
  mainV03_56_B.BusCreator.flightCondition = rtb_FromWs[5];
  mainV03_56_B.BusCreator.Position_cmd[0] =
    mainV03_56_P.ComandedpositioninEarthaxis_Value[0];
  mainV03_56_B.BusCreator.Position_cmd[1] =
    mainV03_56_P.ComandedpositioninEarthaxis_Value[1];
  mainV03_56_B.BusCreator.Position_cmd[2] =
    mainV03_56_P.ComandedpositioninEarthaxis_Value[2];
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Level2 S-Function Block: '<S5>/Joystick Input' (joyinput) */
    {
      SimStruct *rts = mainV03_56_M->childSfunctions[0];
      sfcnOutputs(rts, 1);
    }
  }

  /* BusCreator: '<S5>/Bus Creator1' incorporates:
   *  Constant: '<S5>/Comanded position in Earth axis'
   */
  mainV03_56_B.BusCreator1.roll_cmd = mainV03_56_B.axes[0];
  mainV03_56_B.BusCreator1.pitch_cmd = mainV03_56_B.axes[1];
  mainV03_56_B.BusCreator1.yaw_cmd = mainV03_56_B.axes[2];
  mainV03_56_B.BusCreator1.altitude_cmd = mainV03_56_B.axes[3];
  mainV03_56_B.BusCreator1.Throttle1_cmd = rtb_FromWs[4];
  mainV03_56_B.BusCreator1.flightCondition = rtb_FromWs[5];
  mainV03_56_B.BusCreator1.Position_cmd[0] =
    mainV03_56_P.ComandedpositioninEarthaxis_Value[0];
  mainV03_56_B.BusCreator1.Position_cmd[1] =
    mainV03_56_P.ComandedpositioninEarthaxis_Value[1];
  mainV03_56_B.BusCreator1.Position_cmd[2] =
    mainV03_56_P.ComandedpositioninEarthaxis_Value[2];

  /* ManualSwitch: '<S5>/Manual Switch' */
  if (mainV03_56_P.ManualSwitch_CurrentSetting_b == 1) {
    mainV03_56_B.ManualSwitch = mainV03_56_B.BusCreator;
  } else {
    mainV03_56_B.ManualSwitch = mainV03_56_B.BusCreator1;
  }

  /* End of ManualSwitch: '<S5>/Manual Switch' */

  /* ManualSwitch: '<S3>/Manual Switch' incorporates:
   *  Constant: '<S3>/Constant1'
   */
  if (mainV03_56_P.ManualSwitch_CurrentSetting_oi == 1) {
    mainV03_56_B.ManualSwitch_m = mainV03_56_P.Constant1_Value_e2;
  } else {
    mainV03_56_B.ManualSwitch_m = mainV03_56_B.ManualSwitch.flightCondition;
  }

  /* End of ManualSwitch: '<S3>/Manual Switch' */

  /* DataTypeConversion: '<S167>/Data Type Conversion' */
  riseValLimit = rt_roundd_snf(mainV03_56_B.ManualSwitch_m);
  if (rtIsNaN(riseValLimit) || rtIsInf(riseValLimit)) {
    riseValLimit = 0.0;
  } else {
    riseValLimit = fmod(riseValLimit, 4.294967296E+9);
  }

  mainV03_56_B.DataTypeConversion = riseValLimit < 0.0 ? -(int32_T)(uint32_T)
    -riseValLimit : (int32_T)(uint32_T)riseValLimit;

  /* End of DataTypeConversion: '<S167>/Data Type Conversion' */

  /* RelationalOperator: '<S228>/Compare' incorporates:
   *  Constant: '<S228>/Constant'
   */
  mainV03_56_B.Compare = (mainV03_56_B.DataTypeConversion ==
    mainV03_56_P.CompareToConstant_const_f);

  /* BusCreator: '<S498>/Bus Creator' */
  mainV03_56_B.Sensors.LLAMeas.LatitudeMeas_deg =
    mainV03_56_B.plantData.LLA.Latitude_deg;
  mainV03_56_B.Sensors.LLAMeas.LongitudeMeas_deg =
    mainV03_56_B.plantData.LLA.Longitude_deg;
  mainV03_56_B.Sensors.LLAMeas.AltitudeMeas_m =
    mainV03_56_B.plantData.LLA.Altitude_m;
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Switch: '<S516>/Switch' incorporates:
     *  Abs: '<S516>/Abs'
     *  Bias: '<S516>/Bias'
     *  Bias: '<S516>/Bias1'
     *  Constant: '<S502>/initial_pos'
     *  Constant: '<S516>/Constant2'
     *  Constant: '<S517>/Constant'
     *  Math: '<S516>/Math Function1'
     *  RelationalOperator: '<S517>/Compare'
     */
    if (fabs(mainV03_56_P.LLAtoFlatEarth_LL0[0]) >
        mainV03_56_P.CompareToConstant_const_i) {
      rtb_Switch = rt_modd_snf(mainV03_56_P.LLAtoFlatEarth_LL0[0] +
        mainV03_56_P.Bias_Bias_i4, mainV03_56_P.Constant2_Value_f) +
        mainV03_56_P.Bias1_Bias_p;
    } else {
      rtb_Switch = mainV03_56_P.LLAtoFlatEarth_LL0[0];
    }

    /* End of Switch: '<S516>/Switch' */

    /* Abs: '<S513>/Abs1' */
    rtb_sincos_o2_i_idx_0 = fabs(rtb_Switch);

    /* Switch: '<S513>/Switch' incorporates:
     *  Bias: '<S513>/Bias'
     *  Bias: '<S513>/Bias1'
     *  Constant: '<S504>/Constant'
     *  Constant: '<S504>/Constant1'
     *  Constant: '<S515>/Constant'
     *  Gain: '<S513>/Gain'
     *  Product: '<S513>/Divide1'
     *  RelationalOperator: '<S515>/Compare'
     *  Switch: '<S504>/Switch1'
     */
    if (rtb_sincos_o2_i_idx_0 > mainV03_56_P.CompareToConstant_const_d) {
      /* Signum: '<S513>/Sign1' */
      if (rtb_Switch < 0.0) {
        rtb_Switch = -1.0;
      } else if (rtb_Switch > 0.0) {
        rtb_Switch = 1.0;
      } else {
        if (rtb_Switch == 0.0) {
          rtb_Switch = 0.0;
        }
      }

      /* End of Signum: '<S513>/Sign1' */
      mainV03_56_B.Switch_e = ((rtb_sincos_o2_i_idx_0 +
        mainV03_56_P.Bias_Bias_pz) * mainV03_56_P.Gain_Gain_bu +
        mainV03_56_P.Bias1_Bias_h) * rtb_Switch;
      rtb_Switch = mainV03_56_P.Constant_Value_h;
    } else {
      mainV03_56_B.Switch_e = rtb_Switch;
      rtb_Switch = mainV03_56_P.Constant1_Value_e;
    }

    /* End of Switch: '<S513>/Switch' */

    /* Sum: '<S504>/Sum' incorporates:
     *  Constant: '<S502>/initial_pos'
     */
    rtb_sincos_o2_i_idx_0 = rtb_Switch + mainV03_56_P.LLAtoFlatEarth_LL0[1];

    /* Abs: '<S514>/Abs' */
    rtb_Switch = fabs(rtb_sincos_o2_i_idx_0);

    /* Switch: '<S514>/Switch' incorporates:
     *  Bias: '<S514>/Bias'
     *  Bias: '<S514>/Bias1'
     *  Constant: '<S514>/Constant2'
     *  Constant: '<S518>/Constant'
     *  Math: '<S514>/Math Function1'
     *  RelationalOperator: '<S518>/Compare'
     */
    if (rtb_Switch > mainV03_56_P.CompareToConstant_const_n3) {
      mainV03_56_B.Switch_n = rt_modd_snf(rtb_sincos_o2_i_idx_0 +
        mainV03_56_P.Bias_Bias_m, mainV03_56_P.Constant2_Value_bb) +
        mainV03_56_P.Bias1_Bias_mf;
    } else {
      mainV03_56_B.Switch_n = rtb_sincos_o2_i_idx_0;
    }

    /* End of Switch: '<S514>/Switch' */
  }

  /* Sum: '<S502>/Sum1' */
  mainV03_56_B.Sum1_d[0] = mainV03_56_B.plantData.LLA.Latitude_deg -
    mainV03_56_B.Switch_e;
  mainV03_56_B.Sum1_d[1] = mainV03_56_B.plantData.LLA.Longitude_deg -
    mainV03_56_B.Switch_n;

  /* Abs: '<S510>/Abs' */
  mainV03_56_B.Abs_k = fabs(mainV03_56_B.Sum1_d[0]);

  /* Switch: '<S510>/Switch' incorporates:
   *  Constant: '<S511>/Constant'
   *  RelationalOperator: '<S511>/Compare'
   */
  if (mainV03_56_B.Abs_k > mainV03_56_P.CompareToConstant_const_na) {
    /* Bias: '<S510>/Bias' */
    mainV03_56_B.Bias_a = mainV03_56_B.Sum1_d[0] + mainV03_56_P.Bias_Bias_k;

    /* Bias: '<S510>/Bias1' incorporates:
     *  Constant: '<S510>/Constant2'
     *  Math: '<S510>/Math Function1'
     */
    mainV03_56_B.Bias1_h = rt_modd_snf(mainV03_56_B.Bias_a,
      mainV03_56_P.Constant2_Value_lr) + mainV03_56_P.Bias1_Bias_do;
    mainV03_56_B.Switch_ez = mainV03_56_B.Bias1_h;
  } else {
    mainV03_56_B.Switch_ez = mainV03_56_B.Sum1_d[0];
  }

  /* End of Switch: '<S510>/Switch' */

  /* Abs: '<S507>/Abs1' */
  mainV03_56_B.Abs1_j = fabs(mainV03_56_B.Switch_ez);

  /* Switch: '<S507>/Switch' incorporates:
   *  Constant: '<S503>/Constant'
   *  Constant: '<S503>/Constant1'
   *  Constant: '<S509>/Constant'
   *  RelationalOperator: '<S509>/Compare'
   *  Switch: '<S503>/Switch1'
   */
  if (mainV03_56_B.Abs1_j > mainV03_56_P.CompareToConstant_const_j) {
    /* Bias: '<S507>/Bias' */
    mainV03_56_B.Bias_o = mainV03_56_B.Abs1_j + mainV03_56_P.Bias_Bias_pa;

    /* Gain: '<S507>/Gain' */
    mainV03_56_B.Gain_j = mainV03_56_P.Gain_Gain_cn * mainV03_56_B.Bias_o;

    /* Bias: '<S507>/Bias1' */
    mainV03_56_B.Bias1_hf = mainV03_56_B.Gain_j + mainV03_56_P.Bias1_Bias_m;

    /* Signum: '<S507>/Sign1' */
    if (mainV03_56_B.Switch_ez < 0.0) {
      mainV03_56_B.Sign1_h = -1.0;
    } else if (mainV03_56_B.Switch_ez > 0.0) {
      mainV03_56_B.Sign1_h = 1.0;
    } else if (mainV03_56_B.Switch_ez == 0.0) {
      mainV03_56_B.Sign1_h = 0.0;
    } else {
      mainV03_56_B.Sign1_h = mainV03_56_B.Switch_ez;
    }

    /* End of Signum: '<S507>/Sign1' */

    /* Product: '<S507>/Divide1' */
    mainV03_56_B.Divide1 = mainV03_56_B.Sign1_h * mainV03_56_B.Bias1_hf;
    mainV03_56_B.Switch_p = mainV03_56_B.Divide1;
    mainV03_56_B.Switch1_g = mainV03_56_P.Constant_Value_pu;
  } else {
    mainV03_56_B.Switch_p = mainV03_56_B.Switch_ez;
    mainV03_56_B.Switch1_g = mainV03_56_P.Constant1_Value_g;
  }

  /* End of Switch: '<S507>/Switch' */

  /* Sum: '<S503>/Sum' */
  mainV03_56_B.Sum_oa = mainV03_56_B.Switch1_g + mainV03_56_B.Sum1_d[1];

  /* Abs: '<S508>/Abs' */
  mainV03_56_B.Abs_d = fabs(mainV03_56_B.Sum_oa);

  /* Switch: '<S508>/Switch' incorporates:
   *  Constant: '<S512>/Constant'
   *  RelationalOperator: '<S512>/Compare'
   */
  if (mainV03_56_B.Abs_d > mainV03_56_P.CompareToConstant_const_k) {
    /* Bias: '<S508>/Bias' */
    mainV03_56_B.Bias = mainV03_56_B.Sum_oa + mainV03_56_P.Bias_Bias_ok;

    /* Bias: '<S508>/Bias1' incorporates:
     *  Constant: '<S508>/Constant2'
     *  Math: '<S508>/Math Function1'
     */
    mainV03_56_B.Bias1 = rt_modd_snf(mainV03_56_B.Bias,
      mainV03_56_P.Constant2_Value_cg) + mainV03_56_P.Bias1_Bias_i;
    mainV03_56_B.Switch_aa = mainV03_56_B.Bias1;
  } else {
    mainV03_56_B.Switch_aa = mainV03_56_B.Sum_oa;
  }

  /* End of Switch: '<S508>/Switch' */
  /* Unit Conversion - from: deg to: rad
     Expression: output = (0.0174533*input) + (0) */
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Sum: '<S523>/Sum' incorporates:
     *  Constant: '<S523>/Constant'
     *  Constant: '<S523>/f'
     */
    rtb_Switch = mainV03_56_P.f_Value_n - mainV03_56_P.Constant_Value_gt;

    /* Sqrt: '<S524>/sqrt' incorporates:
     *  Constant: '<S524>/Constant'
     *  Product: '<S524>/Product1'
     *  Sum: '<S524>/Sum1'
     */
    rtb_Switch = sqrt(mainV03_56_P.Constant_Value_bf - rtb_Switch * rtb_Switch);

    /* UnitConversion: '<S521>/Unit Conversion' */
    /* Unit Conversion - from: deg to: rad
       Expression: output = (0.0174533*input) + (0) */
    rtb_Abs1 = 0.017453292519943295 * mainV03_56_B.Switch_e;

    /* Trigonometry: '<S522>/Trigonometric Function1' */
    rtb_Sum1_no = sin(rtb_Abs1);

    /* Sum: '<S522>/Sum1' incorporates:
     *  Constant: '<S522>/Constant'
     *  Product: '<S522>/Product1'
     */
    rtb_Sum1_no = mainV03_56_P.Constant_Value_ke - rtb_Switch * rtb_Switch *
      rtb_Sum1_no * rtb_Sum1_no;

    /* Product: '<S520>/Product1' incorporates:
     *  Constant: '<S520>/Constant1'
     *  Sqrt: '<S520>/sqrt'
     */
    rtb_Rn = mainV03_56_P.Constant1_Value_m / sqrt(rtb_Sum1_no);

    /* Trigonometry: '<S520>/Trigonometric Function1' incorporates:
     *  Constant: '<S520>/Constant'
     *  Constant: '<S520>/Constant2'
     *  Product: '<S520>/Product2'
     *  Product: '<S520>/Product3'
     *  Sum: '<S520>/Sum1'
     */
    mainV03_56_B.TrigonometricFunction1_k = rt_atan2d_snf
      (mainV03_56_P.Constant2_Value_n, (mainV03_56_P.Constant_Value_pj -
        rtb_Switch * rtb_Switch) * rtb_Rn / rtb_Sum1_no);

    /* UnitConversion: '<S519>/Unit Conversion' incorporates:
     *  Constant: '<S505>/ref_pos'
     */
    /* Unit Conversion - from: deg to: rad
       Expression: output = (0.0174533*input) + (0) */
    rtb_Switch = 0.017453292519943295 * mainV03_56_P.LLAtoFlatEarth_psi;

    /* Trigonometry: '<S505>/SinCos' */
    mainV03_56_B.SinCos_o1_d = sin(rtb_Switch);
    mainV03_56_B.SinCos_o2_p = cos(rtb_Switch);
  }

  /* Product: '<S505>/dNorth' incorporates:
   *  UnitConversion: '<S506>/Unit Conversion'
   */
  mainV03_56_B.dNorth = 0.017453292519943295 * mainV03_56_B.Switch_p /
    mainV03_56_B.TrigonometricFunction1_k;

  /* Product: '<S505>/x*cos' */
  mainV03_56_B.xcos_c = mainV03_56_B.dNorth * mainV03_56_B.SinCos_o2_p;
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Trigonometry: '<S520>/Trigonometric Function2' incorporates:
     *  Constant: '<S520>/Constant3'
     *  Product: '<S520>/Product4'
     *  Trigonometry: '<S520>/Trigonometric Function'
     */
    mainV03_56_B.TrigonometricFunction2_l = rt_atan2d_snf
      (mainV03_56_P.Constant3_Value_f, rtb_Rn * cos(rtb_Abs1));
  }

  /* Product: '<S505>/dEast' incorporates:
   *  UnitConversion: '<S506>/Unit Conversion'
   */
  mainV03_56_B.dEast = 1.0 / mainV03_56_B.TrigonometricFunction2_l *
    (0.017453292519943295 * mainV03_56_B.Switch_aa);

  /* Product: '<S505>/y*sin' */
  mainV03_56_B.ysin_l = mainV03_56_B.dEast * mainV03_56_B.SinCos_o1_d;

  /* Sum: '<S505>/Sum2' */
  mainV03_56_B.Px = mainV03_56_B.xcos_c + mainV03_56_B.ysin_l;

  /* Product: '<S505>/x*sin' */
  mainV03_56_B.xsin_a = mainV03_56_B.dNorth * mainV03_56_B.SinCos_o1_d;

  /* Product: '<S505>/y*cos' */
  mainV03_56_B.ycos_b = mainV03_56_B.dEast * mainV03_56_B.SinCos_o2_p;

  /* Sum: '<S505>/Sum3' */
  mainV03_56_B.Py = mainV03_56_B.ycos_b - mainV03_56_B.xsin_a;

  /* Sum: '<S502>/Sum' incorporates:
   *  Constant: '<S498>/Constant13'
   */
  mainV03_56_B.alt_p = mainV03_56_B.plantData.LLA.Altitude_m +
    mainV03_56_P.Constant13_Value_i;

  /* UnaryMinus: '<S502>/Ze2height' */
  rtb_UnaryMinus_j = -mainV03_56_B.alt_p;

  /* SignalConversion: '<S490>/TmpSignal ConversionAtBus Creator2Inport2' */
  mainV03_56_B.X_nedMeas[0] = mainV03_56_B.Px;
  mainV03_56_B.X_nedMeas[1] = mainV03_56_B.Py;
  mainV03_56_B.X_nedMeas[2] = rtb_UnaryMinus_j;

  /* Derivative: '<S498>/Derivative' */
  if ((mainV03_56_DW.TimeStampA_p >= mainV03_56_M->Timing.t[0]) &&
      (mainV03_56_DW.TimeStampB_k >= mainV03_56_M->Timing.t[0])) {
    mainV03_56_B.Sensors.VearthMeas[0] = 0.0;
    mainV03_56_B.Sensors.VearthMeas[1] = 0.0;
    mainV03_56_B.Sensors.VearthMeas[2] = 0.0;
  } else {
    rtb_Sum1_no = mainV03_56_DW.TimeStampA_p;
    lastU_0 = (real_T (*)[3])mainV03_56_DW.LastUAtTimeA_g;
    if (mainV03_56_DW.TimeStampA_p < mainV03_56_DW.TimeStampB_k) {
      if (mainV03_56_DW.TimeStampB_k < mainV03_56_M->Timing.t[0]) {
        rtb_Sum1_no = mainV03_56_DW.TimeStampB_k;
        lastU_0 = (real_T (*)[3])mainV03_56_DW.LastUAtTimeB_a;
      }
    } else {
      if (mainV03_56_DW.TimeStampA_p >= mainV03_56_M->Timing.t[0]) {
        rtb_Sum1_no = mainV03_56_DW.TimeStampB_k;
        lastU_0 = (real_T (*)[3])mainV03_56_DW.LastUAtTimeB_a;
      }
    }

    rtb_sincos_o1_n_idx_1 = mainV03_56_M->Timing.t[0] - rtb_Sum1_no;
    mainV03_56_B.Sensors.VearthMeas[0] = (mainV03_56_B.X_nedMeas[0] - (*lastU_0)
      [0]) / rtb_sincos_o1_n_idx_1;
    mainV03_56_B.Sensors.VearthMeas[1] = (mainV03_56_B.X_nedMeas[1] - (*lastU_0)
      [1]) / rtb_sincos_o1_n_idx_1;
    mainV03_56_B.Sensors.VearthMeas[2] = (mainV03_56_B.X_nedMeas[2] - (*lastU_0)
      [2]) / rtb_sincos_o1_n_idx_1;
  }

  /* End of Derivative: '<S498>/Derivative' */

  /* TransferFcn: '<S538>/Transfer Fcn X' */
  mainV03_56_B.TransferFcnX = 0.0;
  mainV03_56_B.TransferFcnX += mainV03_56_P.TransferFcnX_C[0] *
    mainV03_56_X.TransferFcnX_CSTATE[0];
  mainV03_56_B.TransferFcnX += mainV03_56_P.TransferFcnX_C[1] *
    mainV03_56_X.TransferFcnX_CSTATE[1];

  /* TransferFcn: '<S538>/Transfer Fcn Y' */
  mainV03_56_B.TransferFcnY = 0.0;
  mainV03_56_B.TransferFcnY += mainV03_56_P.TransferFcnY_C[0] *
    mainV03_56_X.TransferFcnY_CSTATE[0];
  mainV03_56_B.TransferFcnY += mainV03_56_P.TransferFcnY_C[1] *
    mainV03_56_X.TransferFcnY_CSTATE[1];

  /* TransferFcn: '<S538>/Transfer Fcn Z' */
  mainV03_56_B.TransferFcnZ = 0.0;
  mainV03_56_B.TransferFcnZ += mainV03_56_P.TransferFcnZ_C[0] *
    mainV03_56_X.TransferFcnZ_CSTATE[0];
  mainV03_56_B.TransferFcnZ += mainV03_56_P.TransferFcnZ_C[1] *
    mainV03_56_X.TransferFcnZ_CSTATE[1];

  /* Gain: '<S531>/Zero-Order Hold' */
  mainV03_56_B.ZeroOrderHold[0] = mainV03_56_P.ZeroOrderHold_Gain *
    mainV03_56_B.plantData.Omega_body.p;
  mainV03_56_B.ZeroOrderHold[1] = mainV03_56_P.ZeroOrderHold_Gain *
    mainV03_56_B.plantData.Omega_body.q;
  mainV03_56_B.ZeroOrderHold[2] = mainV03_56_P.ZeroOrderHold_Gain *
    mainV03_56_B.plantData.Omega_body.r;
  for (s155_iter = 0; s155_iter < 3; s155_iter++) {
    /* Gain: '<S531>/Zero-Order Hold1' */
    mainV03_56_B.ZeroOrderHold1[s155_iter] = mainV03_56_P.ZeroOrderHold1_Gain *
      mainV03_56_B.plantData.Accel_body[s155_iter];

    /* Product: '<S499>/Product' */
    mainV03_56_B.Product_b[s155_iter] = 0.0;
    mainV03_56_B.Product_b[s155_iter] +=
      mainV03_56_B.plantData.DCM_body_earth[s155_iter] *
      mainV03_56_B.envData.Gravity_ned[0];
    mainV03_56_B.Product_b[s155_iter] +=
      mainV03_56_B.plantData.DCM_body_earth[s155_iter + 3] *
      mainV03_56_B.envData.Gravity_ned[1];
    mainV03_56_B.Product_b[s155_iter] +=
      mainV03_56_B.plantData.DCM_body_earth[s155_iter + 6] *
      mainV03_56_B.envData.Gravity_ned[2];

    /* Gain: '<S531>/Zero-Order Hold2' */
    mainV03_56_B.ZeroOrderHold2[s155_iter] = mainV03_56_P.ZeroOrderHold2_Gain *
      mainV03_56_B.Product_b[s155_iter];

    /* Gain: '<S531>/Zero-Order Hold4' */
    mainV03_56_B.ZeroOrderHold4[s155_iter] = mainV03_56_P.ZeroOrderHold4_Gain *
      mainV03_56_B.plantData.CG[s155_iter];

    /* Sum: '<S531>/Sum7' incorporates:
     *  Constant: '<S531>/wl_ins'
     */
    mainV03_56_B.Sum7[s155_iter] = mainV03_56_B.ZeroOrderHold4[s155_iter] -
      mainV03_56_P.ThreeaxisInertialMeasurementUnit_imu[s155_iter];

    /* Gain: '<S531>/Gain' */
    mainV03_56_B.Gain_d[s155_iter] = mainV03_56_P.Gain_Gain_iwz[s155_iter] *
      mainV03_56_B.Sum7[s155_iter];

    /* Gain: '<S531>/Zero-Order Hold3' */
    mainV03_56_B.ZeroOrderHold3[s155_iter] = mainV03_56_P.ZeroOrderHold3_Gain *
      mainV03_56_B.plantData.dOmega_body[s155_iter];
  }

  /* Product: '<S543>/j x k' */
  mainV03_56_B.jxk_b = mainV03_56_B.ZeroOrderHold[1] * mainV03_56_B.Gain_d[2];

  /* Product: '<S543>/k x i' */
  mainV03_56_B.kxi_f = mainV03_56_B.ZeroOrderHold[2] * mainV03_56_B.Gain_d[0];

  /* Product: '<S543>/i x j' */
  mainV03_56_B.ixj_e = mainV03_56_B.ZeroOrderHold[0] * mainV03_56_B.Gain_d[1];

  /* Product: '<S544>/k x j' */
  mainV03_56_B.kxj_l = mainV03_56_B.ZeroOrderHold[2] * mainV03_56_B.Gain_d[1];

  /* Product: '<S544>/i x k' */
  mainV03_56_B.ixk_o2 = mainV03_56_B.ZeroOrderHold[0] * mainV03_56_B.Gain_d[2];

  /* Product: '<S544>/j x i' */
  mainV03_56_B.jxi_fh = mainV03_56_B.ZeroOrderHold[1] * mainV03_56_B.Gain_d[0];

  /* Sum: '<S540>/Sum' */
  mainV03_56_B.Sum_b[0] = mainV03_56_B.jxk_b - mainV03_56_B.kxj_l;
  mainV03_56_B.Sum_b[1] = mainV03_56_B.kxi_f - mainV03_56_B.ixk_o2;
  mainV03_56_B.Sum_b[2] = mainV03_56_B.ixj_e - mainV03_56_B.jxi_fh;

  /* Product: '<S541>/j x k' */
  mainV03_56_B.jxk_cj = mainV03_56_B.ZeroOrderHold[1] * mainV03_56_B.Sum_b[2];

  /* Product: '<S541>/k x i' */
  mainV03_56_B.kxi_l1 = mainV03_56_B.ZeroOrderHold[2] * mainV03_56_B.Sum_b[0];

  /* Product: '<S541>/i x j' */
  mainV03_56_B.ixj_k = mainV03_56_B.ZeroOrderHold[0] * mainV03_56_B.Sum_b[1];

  /* Product: '<S542>/k x j' */
  mainV03_56_B.kxj_me = mainV03_56_B.ZeroOrderHold[2] * mainV03_56_B.Sum_b[1];

  /* Product: '<S542>/i x k' */
  mainV03_56_B.ixk_fv = mainV03_56_B.ZeroOrderHold[0] * mainV03_56_B.Sum_b[2];

  /* Product: '<S542>/j x i' */
  mainV03_56_B.jxi_j3 = mainV03_56_B.ZeroOrderHold[1] * mainV03_56_B.Sum_b[0];

  /* Sum: '<S539>/Sum' */
  mainV03_56_B.Sum_h[0] = mainV03_56_B.jxk_cj - mainV03_56_B.kxj_me;
  mainV03_56_B.Sum_h[1] = mainV03_56_B.kxi_l1 - mainV03_56_B.ixk_fv;
  mainV03_56_B.Sum_h[2] = mainV03_56_B.ixj_k - mainV03_56_B.jxi_j3;

  /* Product: '<S545>/j x k' */
  mainV03_56_B.jxk_l = mainV03_56_B.ZeroOrderHold3[1] * mainV03_56_B.Gain_d[2];

  /* Product: '<S545>/k x i' */
  mainV03_56_B.kxi_b = mainV03_56_B.ZeroOrderHold3[2] * mainV03_56_B.Gain_d[0];

  /* Product: '<S545>/i x j' */
  mainV03_56_B.ixj_d = mainV03_56_B.ZeroOrderHold3[0] * mainV03_56_B.Gain_d[1];

  /* Product: '<S546>/k x j' */
  mainV03_56_B.kxj_ly = mainV03_56_B.ZeroOrderHold3[2] * mainV03_56_B.Gain_d[1];

  /* Product: '<S546>/i x k' */
  mainV03_56_B.ixk_h = mainV03_56_B.ZeroOrderHold3[0] * mainV03_56_B.Gain_d[2];

  /* Product: '<S546>/j x i' */
  mainV03_56_B.jxi_ku = mainV03_56_B.ZeroOrderHold3[1] * mainV03_56_B.Gain_d[0];

  /* Sum: '<S536>/Sum' */
  mainV03_56_B.Sum_m[0] = mainV03_56_B.jxk_l - mainV03_56_B.kxj_ly;
  mainV03_56_B.Sum_m[1] = mainV03_56_B.kxi_b - mainV03_56_B.ixk_h;
  mainV03_56_B.Sum_m[2] = mainV03_56_B.ixj_d - mainV03_56_B.jxi_ku;

  /* Sum: '<S531>/Sum' */
  mainV03_56_B.Sum_cg[0] = ((mainV03_56_B.ZeroOrderHold1[0] -
    mainV03_56_B.ZeroOrderHold2[0]) + mainV03_56_B.Sum_h[0]) +
    mainV03_56_B.Sum_m[0];
  mainV03_56_B.Sum_cg[1] = ((mainV03_56_B.ZeroOrderHold1[1] -
    mainV03_56_B.ZeroOrderHold2[1]) + mainV03_56_B.Sum_h[1]) +
    mainV03_56_B.Sum_m[1];
  mainV03_56_B.Sum_cg[2] = ((mainV03_56_B.ZeroOrderHold1[2] -
    mainV03_56_B.ZeroOrderHold2[2]) + mainV03_56_B.Sum_h[2]) +
    mainV03_56_B.Sum_m[2];
  for (s155_iter = 0; s155_iter < 3; s155_iter++) {
    /* Product: '<S531>/Product' incorporates:
     *  Constant: '<S531>/Scale factors & Cross-coupling  errors'
     */
    mainV03_56_B.Product_g[s155_iter] = 0.0;
    mainV03_56_B.Product_g[s155_iter] +=
      mainV03_56_P.ThreeaxisInertialMeasurementUnit_a_sf_cc[s155_iter] *
      mainV03_56_B.Sum_cg[0];
    mainV03_56_B.Product_g[s155_iter] +=
      mainV03_56_P.ThreeaxisInertialMeasurementUnit_a_sf_cc[s155_iter + 3] *
      mainV03_56_B.Sum_cg[1];
    mainV03_56_B.Product_g[s155_iter] +=
      mainV03_56_P.ThreeaxisInertialMeasurementUnit_a_sf_cc[s155_iter + 6] *
      mainV03_56_B.Sum_cg[2];

    /* Sum: '<S531>/Sum4' incorporates:
     *  Constant: '<S531>/Measurement bias'
     */
    mainV03_56_B.Sum4_d[s155_iter] = mainV03_56_B.Product_g[s155_iter] +
      mainV03_56_P.ThreeaxisInertialMeasurementUnit_a_bias[s155_iter];
  }

  /* Switch: '<S533>/Switch' incorporates:
   *  Constant: '<S533>/Constant'
   */
  if (mainV03_56_P.Constant_Value_ob >= mainV03_56_P.Switch_Threshold_o5) {
    mainV03_56_B.Switch_o[0] = mainV03_56_B.TransferFcnX;
    mainV03_56_B.Switch_o[1] = mainV03_56_B.TransferFcnY;
    mainV03_56_B.Switch_o[2] = mainV03_56_B.TransferFcnZ;
  } else {
    mainV03_56_B.Switch_o[0] = mainV03_56_B.Sum4_d[0];
    mainV03_56_B.Switch_o[1] = mainV03_56_B.Sum4_d[1];
    mainV03_56_B.Switch_o[2] = mainV03_56_B.Sum4_d[2];
  }

  /* End of Switch: '<S533>/Switch' */
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[2] == 0) {
    /* Gain: '<S534>/Output' incorporates:
     *  RandomNumber: '<S534>/White Noise'
     */
    mainV03_56_B.Output_e[0] = mainV03_56_P.Output_Gain_o[0] *
      mainV03_56_DW.NextOutput_e[0];
    mainV03_56_B.Output_e[1] = mainV03_56_P.Output_Gain_o[1] *
      mainV03_56_DW.NextOutput_e[1];
    mainV03_56_B.Output_e[2] = mainV03_56_P.Output_Gain_o[2] *
      mainV03_56_DW.NextOutput_e[2];
  }

  /* Sum: '<S531>/Sum1' */
  mainV03_56_B.Sum1_e[0] = mainV03_56_B.Switch_o[0] + mainV03_56_B.Output_e[0];

  /* Saturate: '<S531>/Saturation' */
  if (mainV03_56_B.Sum1_e[0] > mainV03_56_P.Saturation_UpperSat_a[0]) {
    mainV03_56_B.Sensors.AccelMeas_body[0] = mainV03_56_P.Saturation_UpperSat_a
      [0];
  } else if (mainV03_56_B.Sum1_e[0] < mainV03_56_P.Saturation_LowerSat_o[0]) {
    mainV03_56_B.Sensors.AccelMeas_body[0] = mainV03_56_P.Saturation_LowerSat_o
      [0];
  } else {
    mainV03_56_B.Sensors.AccelMeas_body[0] = mainV03_56_B.Sum1_e[0];
  }

  /* Sum: '<S531>/Sum1' */
  mainV03_56_B.Sum1_e[1] = mainV03_56_B.Switch_o[1] + mainV03_56_B.Output_e[1];

  /* Saturate: '<S531>/Saturation' */
  if (mainV03_56_B.Sum1_e[1] > mainV03_56_P.Saturation_UpperSat_a[1]) {
    mainV03_56_B.Sensors.AccelMeas_body[1] = mainV03_56_P.Saturation_UpperSat_a
      [1];
  } else if (mainV03_56_B.Sum1_e[1] < mainV03_56_P.Saturation_LowerSat_o[1]) {
    mainV03_56_B.Sensors.AccelMeas_body[1] = mainV03_56_P.Saturation_LowerSat_o
      [1];
  } else {
    mainV03_56_B.Sensors.AccelMeas_body[1] = mainV03_56_B.Sum1_e[1];
  }

  /* Sum: '<S531>/Sum1' */
  mainV03_56_B.Sum1_e[2] = mainV03_56_B.Switch_o[2] + mainV03_56_B.Output_e[2];

  /* Saturate: '<S531>/Saturation' */
  if (mainV03_56_B.Sum1_e[2] > mainV03_56_P.Saturation_UpperSat_a[2]) {
    mainV03_56_B.Sensors.AccelMeas_body[2] = mainV03_56_P.Saturation_UpperSat_a
      [2];
  } else if (mainV03_56_B.Sum1_e[2] < mainV03_56_P.Saturation_LowerSat_o[2]) {
    mainV03_56_B.Sensors.AccelMeas_body[2] = mainV03_56_P.Saturation_LowerSat_o
      [2];
  } else {
    mainV03_56_B.Sensors.AccelMeas_body[2] = mainV03_56_B.Sum1_e[2];
  }

  /* TransferFcn: '<S550>/Transfer Fcn X' */
  mainV03_56_B.TransferFcnX_k = 0.0;
  mainV03_56_B.TransferFcnX_k += mainV03_56_P.TransferFcnX_C_d[0] *
    mainV03_56_X.TransferFcnX_CSTATE_e[0];
  mainV03_56_B.TransferFcnX_k += mainV03_56_P.TransferFcnX_C_d[1] *
    mainV03_56_X.TransferFcnX_CSTATE_e[1];

  /* TransferFcn: '<S550>/Transfer Fcn Y' */
  mainV03_56_B.TransferFcnY_c = 0.0;
  mainV03_56_B.TransferFcnY_c += mainV03_56_P.TransferFcnY_C_p[0] *
    mainV03_56_X.TransferFcnY_CSTATE_j[0];
  mainV03_56_B.TransferFcnY_c += mainV03_56_P.TransferFcnY_C_p[1] *
    mainV03_56_X.TransferFcnY_CSTATE_j[1];

  /* TransferFcn: '<S550>/Transfer Fcn Z' */
  mainV03_56_B.TransferFcnZ_b = 0.0;
  mainV03_56_B.TransferFcnZ_b += mainV03_56_P.TransferFcnZ_C_d[0] *
    mainV03_56_X.TransferFcnZ_CSTATE_l[0];
  mainV03_56_B.TransferFcnZ_b += mainV03_56_P.TransferFcnZ_C_d[1] *
    mainV03_56_X.TransferFcnZ_CSTATE_l[1];

  /* Gain: '<S532>/Zero-Order Hold' */
  mainV03_56_B.ZeroOrderHold_m[0] = mainV03_56_P.ZeroOrderHold_Gain_c *
    mainV03_56_B.plantData.Omega_body.p;
  mainV03_56_B.ZeroOrderHold_m[1] = mainV03_56_P.ZeroOrderHold_Gain_c *
    mainV03_56_B.plantData.Omega_body.q;
  mainV03_56_B.ZeroOrderHold_m[2] = mainV03_56_P.ZeroOrderHold_Gain_c *
    mainV03_56_B.plantData.Omega_body.r;

  /* Unit Conversion - from: m/s^2 to: gn
     Expression: output = (0.101972*input) + (0) */
  for (s155_iter = 0; s155_iter < 3; s155_iter++) {
    /* Product: '<S532>/Product' incorporates:
     *  Constant: '<S532>/Scale factors & Cross-coupling  errors '
     */
    mainV03_56_B.Product_a[s155_iter] = 0.0;
    mainV03_56_B.Product_a[s155_iter] +=
      mainV03_56_P.ThreeaxisInertialMeasurementUnit_g_sf_cc[s155_iter] *
      mainV03_56_B.ZeroOrderHold_m[0];
    mainV03_56_B.Product_a[s155_iter] +=
      mainV03_56_P.ThreeaxisInertialMeasurementUnit_g_sf_cc[s155_iter + 3] *
      mainV03_56_B.ZeroOrderHold_m[1];
    mainV03_56_B.Product_a[s155_iter] +=
      mainV03_56_P.ThreeaxisInertialMeasurementUnit_g_sf_cc[s155_iter + 6] *
      mainV03_56_B.ZeroOrderHold_m[2];

    /* Gain: '<S532>/Zero-Order Hold1' incorporates:
     *  UnitConversion: '<S530>/Unit Conversion'
     */
    mainV03_56_B.ZeroOrderHold1_f[s155_iter] = 0.10197162129779282 *
      mainV03_56_B.plantData.Accel_body[s155_iter] *
      mainV03_56_P.ZeroOrderHold1_Gain_k;

    /* Product: '<S532>/Product1' incorporates:
     *  Constant: '<S532>/g-sensitive bias'
     */
    mainV03_56_B.Product1_aj[s155_iter] =
      mainV03_56_B.ZeroOrderHold1_f[s155_iter] *
      mainV03_56_P.ThreeaxisInertialMeasurementUnit_g_sens[s155_iter];

    /* Sum: '<S532>/Sum4' incorporates:
     *  Constant: '<S532>/Measurement bias'
     */
    mainV03_56_B.Sum4_a[s155_iter] = (mainV03_56_B.Product_a[s155_iter] +
      mainV03_56_P.ThreeaxisInertialMeasurementUnit_g_bias[s155_iter]) +
      mainV03_56_B.Product1_aj[s155_iter];
  }

  /* Switch: '<S547>/Switch' incorporates:
   *  Constant: '<S547>/Constant'
   */
  if (mainV03_56_P.Constant_Value_id >= mainV03_56_P.Switch_Threshold_f) {
    mainV03_56_B.Switch_l[0] = mainV03_56_B.TransferFcnX_k;
    mainV03_56_B.Switch_l[1] = mainV03_56_B.TransferFcnY_c;
    mainV03_56_B.Switch_l[2] = mainV03_56_B.TransferFcnZ_b;
  } else {
    mainV03_56_B.Switch_l[0] = mainV03_56_B.Sum4_a[0];
    mainV03_56_B.Switch_l[1] = mainV03_56_B.Sum4_a[1];
    mainV03_56_B.Switch_l[2] = mainV03_56_B.Sum4_a[2];
  }

  /* End of Switch: '<S547>/Switch' */
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[2] == 0) {
    /* Gain: '<S548>/Output' incorporates:
     *  RandomNumber: '<S548>/White Noise'
     */
    mainV03_56_B.Output_m[0] = mainV03_56_P.Output_Gain_h[0] *
      mainV03_56_DW.NextOutput_b[0];
    mainV03_56_B.Output_m[1] = mainV03_56_P.Output_Gain_h[1] *
      mainV03_56_DW.NextOutput_b[1];
    mainV03_56_B.Output_m[2] = mainV03_56_P.Output_Gain_h[2] *
      mainV03_56_DW.NextOutput_b[2];
  }

  /* Sum: '<S532>/Sum1' */
  mainV03_56_B.Sum1_i[0] = mainV03_56_B.Switch_l[0] + mainV03_56_B.Output_m[0];

  /* Saturate: '<S532>/Saturation' */
  if (mainV03_56_B.Sum1_i[0] > mainV03_56_P.Saturation_UpperSat_c[0]) {
    mainV03_56_B.Saturation_a[0] = mainV03_56_P.Saturation_UpperSat_c[0];
  } else if (mainV03_56_B.Sum1_i[0] < mainV03_56_P.Saturation_LowerSat_m[0]) {
    mainV03_56_B.Saturation_a[0] = mainV03_56_P.Saturation_LowerSat_m[0];
  } else {
    mainV03_56_B.Saturation_a[0] = mainV03_56_B.Sum1_i[0];
  }

  /* Integrator: '<S525>/phi theta psi' */
  mainV03_56_B.phithetapsi_m[0] = mainV03_56_X.phithetapsi_CSTATE_e[0];

  /* Sum: '<S532>/Sum1' */
  mainV03_56_B.Sum1_i[1] = mainV03_56_B.Switch_l[1] + mainV03_56_B.Output_m[1];

  /* Saturate: '<S532>/Saturation' */
  if (mainV03_56_B.Sum1_i[1] > mainV03_56_P.Saturation_UpperSat_c[1]) {
    mainV03_56_B.Saturation_a[1] = mainV03_56_P.Saturation_UpperSat_c[1];
  } else if (mainV03_56_B.Sum1_i[1] < mainV03_56_P.Saturation_LowerSat_m[1]) {
    mainV03_56_B.Saturation_a[1] = mainV03_56_P.Saturation_LowerSat_m[1];
  } else {
    mainV03_56_B.Saturation_a[1] = mainV03_56_B.Sum1_i[1];
  }

  /* Integrator: '<S525>/phi theta psi' */
  mainV03_56_B.phithetapsi_m[1] = mainV03_56_X.phithetapsi_CSTATE_e[1];

  /* Sum: '<S532>/Sum1' */
  mainV03_56_B.Sum1_i[2] = mainV03_56_B.Switch_l[2] + mainV03_56_B.Output_m[2];

  /* Saturate: '<S532>/Saturation' */
  if (mainV03_56_B.Sum1_i[2] > mainV03_56_P.Saturation_UpperSat_c[2]) {
    mainV03_56_B.Saturation_a[2] = mainV03_56_P.Saturation_UpperSat_c[2];
  } else if (mainV03_56_B.Sum1_i[2] < mainV03_56_P.Saturation_LowerSat_m[2]) {
    mainV03_56_B.Saturation_a[2] = mainV03_56_P.Saturation_LowerSat_m[2];
  } else {
    mainV03_56_B.Saturation_a[2] = mainV03_56_B.Sum1_i[2];
  }

  /* Integrator: '<S525>/phi theta psi' */
  mainV03_56_B.phithetapsi_m[2] = mainV03_56_X.phithetapsi_CSTATE_e[2];

  /* BusCreator: '<S499>/Bus Creator' */
  mainV03_56_B.Sensors.OmegaMeas_body.pMeas = mainV03_56_B.Saturation_a[0];
  mainV03_56_B.Sensors.OmegaMeas_body.qMeas = mainV03_56_B.Saturation_a[1];
  mainV03_56_B.Sensors.OmegaMeas_body.rMeas = mainV03_56_B.Saturation_a[2];

  /* SignalConversion: '<S527>/TmpSignal ConversionAtsincosInport1' */
  mainV03_56_B.TmpSignalConversionAtsincosInport1_p[0] =
    mainV03_56_B.phithetapsi_m[2];
  mainV03_56_B.TmpSignalConversionAtsincosInport1_p[1] =
    mainV03_56_B.phithetapsi_m[1];
  mainV03_56_B.TmpSignalConversionAtsincosInport1_p[2] =
    mainV03_56_B.phithetapsi_m[0];

  /* Trigonometry: '<S527>/sincos' */
  rtb_Shiftright[0] = sin(mainV03_56_B.TmpSignalConversionAtsincosInport1_p[0]);
  rtb_sincos_o2_p[0] = cos(mainV03_56_B.TmpSignalConversionAtsincosInport1_p[0]);

  /* Sum: '<S500>/Sum' */
  mainV03_56_B.Sensors.VaeroMeas_body[0] = mainV03_56_B.plantData.V_body[0] -
    mainV03_56_B.envData.windVelocity[0];

  /* BusCreator: '<S490>/Bus Creator2' */
  mainV03_56_B.Sensors.X_nedMeas[0] = mainV03_56_B.X_nedMeas[0];

  /* Trigonometry: '<S527>/sincos' */
  rtb_Shiftright[1] = sin(mainV03_56_B.TmpSignalConversionAtsincosInport1_p[1]);
  rtb_sincos_o2_p[1] = cos(mainV03_56_B.TmpSignalConversionAtsincosInport1_p[1]);

  /* Sum: '<S500>/Sum' */
  mainV03_56_B.Sensors.VaeroMeas_body[1] = mainV03_56_B.plantData.V_body[1] -
    mainV03_56_B.envData.windVelocity[1];

  /* BusCreator: '<S490>/Bus Creator2' */
  mainV03_56_B.Sensors.X_nedMeas[1] = mainV03_56_B.X_nedMeas[1];

  /* Trigonometry: '<S527>/sincos' */
  rtb_Shiftright[2] = sin(mainV03_56_B.TmpSignalConversionAtsincosInport1_p[2]);
  rtb_sincos_o2_p[2] = cos(mainV03_56_B.TmpSignalConversionAtsincosInport1_p[2]);

  /* Sum: '<S500>/Sum' */
  mainV03_56_B.Sensors.VaeroMeas_body[2] = mainV03_56_B.plantData.V_body[2] -
    mainV03_56_B.envData.windVelocity[2];

  /* BusCreator: '<S490>/Bus Creator2' */
  mainV03_56_B.Sensors.X_nedMeas[2] = mainV03_56_B.X_nedMeas[2];

  /* Fcn: '<S527>/Fcn11' */
  mainV03_56_B.VectorConcatenate_h[0] = rtb_sincos_o2_p[1] * rtb_sincos_o2_p[0];

  /* Fcn: '<S527>/Fcn21' */
  mainV03_56_B.VectorConcatenate_h[1] = rtb_Shiftright[2] * rtb_Shiftright[1] *
    rtb_sincos_o2_p[0] - rtb_sincos_o2_p[2] * rtb_Shiftright[0];

  /* Fcn: '<S527>/Fcn31' */
  mainV03_56_B.VectorConcatenate_h[2] = rtb_sincos_o2_p[2] * rtb_Shiftright[1] *
    rtb_sincos_o2_p[0] + rtb_Shiftright[2] * rtb_Shiftright[0];

  /* Fcn: '<S527>/Fcn12' */
  mainV03_56_B.VectorConcatenate_h[3] = rtb_sincos_o2_p[1] * rtb_Shiftright[0];

  /* Fcn: '<S527>/Fcn22' */
  mainV03_56_B.VectorConcatenate_h[4] = rtb_Shiftright[2] * rtb_Shiftright[1] *
    rtb_Shiftright[0] + rtb_sincos_o2_p[2] * rtb_sincos_o2_p[0];

  /* Fcn: '<S527>/Fcn32' */
  mainV03_56_B.VectorConcatenate_h[5] = rtb_sincos_o2_p[2] * rtb_Shiftright[1] *
    rtb_Shiftright[0] - rtb_Shiftright[2] * rtb_sincos_o2_p[0];

  /* Fcn: '<S527>/Fcn13' */
  mainV03_56_B.VectorConcatenate_h[6] = -rtb_Shiftright[1];

  /* Fcn: '<S527>/Fcn23' */
  mainV03_56_B.VectorConcatenate_h[7] = rtb_Shiftright[2] * rtb_sincos_o2_p[1];

  /* Fcn: '<S527>/Fcn33' */
  mainV03_56_B.VectorConcatenate_h[8] = rtb_sincos_o2_p[2] * rtb_sincos_o2_p[1];

  /* BusCreator: '<S490>/Bus Creator2' */
  memcpy(&mainV03_56_B.Sensors.DCMMeas_body_earth[0],
         &mainV03_56_B.VectorConcatenate_h[0], 9U * sizeof(real_T));
  mainV03_56_B.Sensors.alphaMeas = mainV03_56_B.plantData.alpha;
  mainV03_56_B.Sensors.betaMeas = mainV03_56_B.plantData.beta;
  mainV03_56_B.Sensors.EulerMeas[0] = mainV03_56_B.phithetapsi_m[0];
  mainV03_56_B.Sensors.CGMeas[0] = mainV03_56_B.plantData.CG[0];
  mainV03_56_B.Sensors.EulerMeas[1] = mainV03_56_B.phithetapsi_m[1];
  mainV03_56_B.Sensors.CGMeas[1] = mainV03_56_B.plantData.CG[1];
  mainV03_56_B.Sensors.EulerMeas[2] = mainV03_56_B.phithetapsi_m[2];
  mainV03_56_B.Sensors.CGMeas[2] = mainV03_56_B.plantData.CG[2];
  mainV03_56_B.Sensors.remainingCapacityMeas =
    mainV03_56_B.plantData.remainingCapacity;
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* SignalConversion: '<S3>/HiddenBuf_InsertedFor_Quadcopter_at_inport_2' */
    mainV03_56_B.HiddenBuf_InsertedFor_Quadcopter_at_inport_2 =
      mainV03_56_B.Compare;

    /* Outputs for Enabled SubSystem: '<S3>/Quadcopter' incorporates:
     *  EnablePort: '<S165>/Enable'
     */
    if (rtmIsMajorTimeStep(mainV03_56_M)) {
      if (mainV03_56_B.HiddenBuf_InsertedFor_Quadcopter_at_inport_2) {
        if (!mainV03_56_DW.Quadcopter_MODE) {
          mainV03_56_DW.Quadcopter_MODE = true;
        }
      } else {
        if (mainV03_56_DW.Quadcopter_MODE) {
          mainV03_56_DW.Quadcopter_MODE = false;
        }
      }
    }

    /* End of Outputs for SubSystem: '<S3>/Quadcopter' */
  }

  /* Outputs for Enabled SubSystem: '<S3>/Quadcopter' incorporates:
   *  EnablePort: '<S165>/Enable'
   */
  if (mainV03_56_DW.Quadcopter_MODE) {
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* BusCreator: '<S165>/Bus Creator' incorporates:
       *  Constant: '<S165>/Constant1'
       *  Constant: '<S165>/Constant10'
       *  Constant: '<S165>/Constant2'
       *  Constant: '<S165>/Constant9'
       *  Gain: '<S165>/Gain'
       *  Gain: '<S165>/Gain1'
       *  Gain: '<S165>/Gain2'
       *  Gain: '<S165>/Gain3'
       */
      mainV03_56_B.actuators_h.deltae = mainV03_56_P.Gain_Gain_dd *
        mainV03_56_P.deltae_degrees;
      mainV03_56_B.actuators_h.deltar = mainV03_56_P.Gain1_Gain_e *
        mainV03_56_P.deltar_degrees;
      mainV03_56_B.actuators_h.deltafr = mainV03_56_P.Gain2_Gain_mx *
        mainV03_56_P.deltafr_degrees;
      mainV03_56_B.actuators_h.deltafl = mainV03_56_P.Gain3_Gain_e *
        mainV03_56_P.deltafl_degrees;

      /* BusCreator: '<S165>/Bus Creator1' incorporates:
       *  Constant: '<S165>/Constant20'
       *  Constant: '<S165>/Constant21'
       *  Constant: '<S165>/Constant22'
       *  Constant: '<S165>/Constant23'
       *  Constant: '<S165>/Constant24'
       */
      rtb_sincos_o1_idx_0 = mainV03_56_P.Constant21_Value_f;
      rtb_sincos_o2_g_idx_0 = mainV03_56_P.Constant22_Value_i;
      rtb_sincos_o1_idx_1 = mainV03_56_P.Constant23_Value_p;
      rtb_sincos_o2_g_idx_1 = mainV03_56_P.Constant24_Value_j;
      rtb_sincos_o1_idx_2 = mainV03_56_P.Constant20_Value_j;

      /* BusCreator: '<S165>/Bus Creator3' incorporates:
       *  Constant: '<S165>/deltat_10'
       *  Constant: '<S165>/deltat_6'
       *  Constant: '<S165>/deltat_7'
       *  Constant: '<S165>/deltat_8'
       *  Constant: '<S165>/deltat_9'
       */
      rtb_sincos_o2_g_idx_2 = mainV03_56_P.deltat_1;
      rtb_BusCreator3_n_Throttle2 = mainV03_56_P.deltat_2;
      rtb_BusCreator3_n_Throttle3 = mainV03_56_P.deltat_3;
      rtb_BusCreator3_n_Throttle4 = mainV03_56_P.deltat_4;
      rtb_BusCreator3_n_Throttle5 = mainV03_56_P.deltat_5;
    }

    /* SignalConversion: '<S194>/TmpSignal ConversionAtsincosInport1' incorporates:
     *  Constant: '<S165>/Constant6'
     */
    mainV03_56_B.TmpSignalConversionAtsincosInport1_i[0] =
      mainV03_56_B.Sensors.EulerMeas[2];
    mainV03_56_B.TmpSignalConversionAtsincosInport1_i[1] =
      mainV03_56_P.Constant6_Value_m[0];
    mainV03_56_B.TmpSignalConversionAtsincosInport1_i[2] =
      mainV03_56_P.Constant6_Value_m[1];

    /* Trigonometry: '<S194>/sincos' */
    rtb_Rn = sin(mainV03_56_B.TmpSignalConversionAtsincosInport1_i[0]);
    rtb_Abs1 = cos(mainV03_56_B.TmpSignalConversionAtsincosInport1_i[0]);
    rtb_Sum1_no = sin(mainV03_56_B.TmpSignalConversionAtsincosInport1_i[1]);
    rtb_sincos_o2_i_idx_0 = cos
      (mainV03_56_B.TmpSignalConversionAtsincosInport1_i[1]);
    rtb_sincos_o2_i_idx_1 = sin
      (mainV03_56_B.TmpSignalConversionAtsincosInport1_i[2]);
    rtb_sincos_o1_n_idx_0 = cos
      (mainV03_56_B.TmpSignalConversionAtsincosInport1_i[2]);

    /* Fcn: '<S194>/Fcn11' */
    mainV03_56_B.VectorConcatenate_lt[0] = rtb_sincos_o2_i_idx_0 * rtb_Abs1;

    /* Fcn: '<S194>/Fcn21' */
    mainV03_56_B.VectorConcatenate_lt[1] = rtb_sincos_o2_i_idx_1 * rtb_Sum1_no *
      rtb_Abs1 - rtb_sincos_o1_n_idx_0 * rtb_Rn;

    /* Fcn: '<S194>/Fcn31' */
    mainV03_56_B.VectorConcatenate_lt[2] = rtb_sincos_o1_n_idx_0 * rtb_Sum1_no *
      rtb_Abs1 + rtb_sincos_o2_i_idx_1 * rtb_Rn;

    /* Fcn: '<S194>/Fcn12' */
    mainV03_56_B.VectorConcatenate_lt[3] = rtb_sincos_o2_i_idx_0 * rtb_Rn;

    /* Fcn: '<S194>/Fcn22' */
    mainV03_56_B.VectorConcatenate_lt[4] = rtb_sincos_o2_i_idx_1 * rtb_Sum1_no *
      rtb_Rn + rtb_sincos_o1_n_idx_0 * rtb_Abs1;

    /* Fcn: '<S194>/Fcn32' */
    mainV03_56_B.VectorConcatenate_lt[5] = rtb_sincos_o1_n_idx_0 * rtb_Sum1_no *
      rtb_Rn - rtb_sincos_o2_i_idx_1 * rtb_Abs1;

    /* Fcn: '<S194>/Fcn13' */
    mainV03_56_B.VectorConcatenate_lt[6] = -rtb_Sum1_no;

    /* Fcn: '<S194>/Fcn23' */
    mainV03_56_B.VectorConcatenate_lt[7] = rtb_sincos_o2_i_idx_1 *
      rtb_sincos_o2_i_idx_0;

    /* Fcn: '<S194>/Fcn33' */
    mainV03_56_B.VectorConcatenate_lt[8] = rtb_sincos_o1_n_idx_0 *
      rtb_sincos_o2_i_idx_0;
    for (rtb_Sum1_eh = 0; rtb_Sum1_eh < 3; rtb_Sum1_eh++) {
      /* Product: '<S165>/To Body Axes3' */
      mainV03_56_B.ToBodyAxes3_l[rtb_Sum1_eh] = 0.0;

      /* Product: '<S165>/To Body Axes4' */
      mainV03_56_B.ToBodyAxes4_p[rtb_Sum1_eh] = 0.0;

      /* Product: '<S165>/To Body Axes3' */
      mainV03_56_B.ToBodyAxes3_l[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_lt[rtb_Sum1_eh] *
        mainV03_56_B.ManualSwitch.Position_cmd[0];

      /* Product: '<S165>/To Body Axes4' */
      mainV03_56_B.ToBodyAxes4_p[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_lt[rtb_Sum1_eh] *
        mainV03_56_B.Sensors.X_nedMeas[0];

      /* Product: '<S165>/To Body Axes3' */
      mainV03_56_B.ToBodyAxes3_l[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_lt[rtb_Sum1_eh + 3] *
        mainV03_56_B.ManualSwitch.Position_cmd[1];

      /* Product: '<S165>/To Body Axes4' */
      mainV03_56_B.ToBodyAxes4_p[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_lt[rtb_Sum1_eh + 3] *
        mainV03_56_B.Sensors.X_nedMeas[1];

      /* Product: '<S165>/To Body Axes3' */
      mainV03_56_B.ToBodyAxes3_l[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_lt[rtb_Sum1_eh + 6] *
        mainV03_56_B.ManualSwitch.Position_cmd[2];

      /* Product: '<S165>/To Body Axes4' */
      mainV03_56_B.ToBodyAxes4_p[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_lt[rtb_Sum1_eh + 6] *
        mainV03_56_B.Sensors.X_nedMeas[2];
    }

    /* ManualSwitch: '<S165>/Manual Switch12' incorporates:
     *  Constant: '<S165>/Constant13'
     */
    if (mainV03_56_P.ManualSwitch12_CurrentSetting == 1) {
      mainV03_56_B.ManualSwitch12_b = mainV03_56_B.ToBodyAxes3_l[0];
    } else {
      mainV03_56_B.ManualSwitch12_b = mainV03_56_P.Constant13_Value_f;
    }

    /* End of ManualSwitch: '<S165>/Manual Switch12' */

    /* Sum: '<S165>/Sum14' */
    mainV03_56_B.Sum14_o = mainV03_56_B.ManualSwitch12_b -
      mainV03_56_B.ToBodyAxes4_p[0];

    /* Gain: '<S179>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_mg = mainV03_56_P.ESTE2_D * mainV03_56_B.Sum14_o;

    /* Integrator: '<S179>/Filter' */
    mainV03_56_B.Filter_et = mainV03_56_X.Filter_CSTATE_ci;

    /* Sum: '<S179>/SumD' */
    mainV03_56_B.SumD_i0 = mainV03_56_B.DerivativeGain_mg -
      mainV03_56_B.Filter_et;

    /* Gain: '<S179>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_ct = mainV03_56_P.ESTE2_N *
      mainV03_56_B.SumD_i0;

    /* Gain: '<S179>/Integral Gain' */
    mainV03_56_B.IntegralGain_mk = mainV03_56_P.ESTE2_I * mainV03_56_B.Sum14_o;

    /* Integrator: '<S179>/Integrator' */
    mainV03_56_B.Integrator_d5 = mainV03_56_X.Integrator_CSTATE_nu;

    /* Gain: '<S179>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_hd = mainV03_56_P.ESTE2_P *
      mainV03_56_B.Sum14_o;

    /* Sum: '<S179>/Sum' */
    mainV03_56_B.Sum_pnb = (mainV03_56_B.ProportionalGain_hd +
      mainV03_56_B.Integrator_d5) + mainV03_56_B.FilterCoefficient_ct;

    /* Sum: '<S165>/Sum2' */
    mainV03_56_B.delta_p = mainV03_56_B.ToBodyAxes3_l[2] -
      mainV03_56_B.Sensors.LLAMeas.AltitudeMeas_m;

    /* Saturate: '<S165>/Saturation9' */
    if (mainV03_56_B.delta_p > mainV03_56_P.Saturation9_UpperSat) {
      mainV03_56_B.Saturation9_k = mainV03_56_P.Saturation9_UpperSat;
    } else if (mainV03_56_B.delta_p < mainV03_56_P.Saturation9_LowerSat) {
      mainV03_56_B.Saturation9_k = mainV03_56_P.Saturation9_LowerSat;
    } else {
      mainV03_56_B.Saturation9_k = mainV03_56_B.delta_p;
    }

    /* End of Saturate: '<S165>/Saturation9' */

    /* Gain: '<S184>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_dp = mainV03_56_P.PIDHeightController1_P *
      mainV03_56_B.Saturation9_k;

    /* Integrator: '<S184>/Integrator' */
    mainV03_56_B.Integrator_n4 = mainV03_56_X.Integrator_CSTATE_aw;

    /* Gain: '<S184>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_j4 = mainV03_56_P.PIDHeightController1_D *
      mainV03_56_B.Saturation9_k;

    /* Integrator: '<S184>/Filter' */
    mainV03_56_B.Filter_l3 = mainV03_56_X.Filter_CSTATE_aw;

    /* Sum: '<S184>/SumD' */
    mainV03_56_B.SumD_j = mainV03_56_B.DerivativeGain_j4 -
      mainV03_56_B.Filter_l3;

    /* Gain: '<S184>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_dn = mainV03_56_P.PIDHeightController1_N *
      mainV03_56_B.SumD_j;
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* Saturate: '<S165>/Saturation3' incorporates:
       *  Constant: '<S165>/Constant'
       */
      if (mainV03_56_P.Constant_Value_d > mainV03_56_P.Saturation3_UpperSat) {
        mainV03_56_B.Saturation3_o = mainV03_56_P.Saturation3_UpperSat;
      } else if (mainV03_56_P.Constant_Value_d <
                 mainV03_56_P.Saturation3_LowerSat) {
        mainV03_56_B.Saturation3_o = mainV03_56_P.Saturation3_LowerSat;
      } else {
        mainV03_56_B.Saturation3_o = mainV03_56_P.Constant_Value_d;
      }

      /* End of Saturate: '<S165>/Saturation3' */
    }

    /* Sum: '<S165>/Sum11' incorporates:
     *  UnaryMinus: '<S165>/Unary Minus1'
     */
    mainV03_56_B.Sum11_j = mainV03_56_B.Saturation3_o -
      (-mainV03_56_B.Sensors.AccelMeas_body[2]);

    /* Gain: '<S180>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_lh = mainV03_56_P.PIDAltitudeAcceleration_P *
      mainV03_56_B.Sum11_j;

    /* Integrator: '<S180>/Integrator' */
    mainV03_56_B.Integrator_bx = mainV03_56_X.Integrator_CSTATE_ij;

    /* Gain: '<S180>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_in = mainV03_56_P.PIDAltitudeAcceleration_D *
      mainV03_56_B.Sum11_j;

    /* Integrator: '<S180>/Filter' */
    mainV03_56_B.Filter_fe = mainV03_56_X.Filter_CSTATE_jg;

    /* Sum: '<S180>/SumD' */
    mainV03_56_B.SumD_ji = mainV03_56_B.DerivativeGain_in -
      mainV03_56_B.Filter_fe;

    /* Gain: '<S180>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_im = mainV03_56_P.PIDAltitudeAcceleration_N *
      mainV03_56_B.SumD_ji;

    /* ManualSwitch: '<S165>/Manual Switch' */
    if (mainV03_56_P.ManualSwitch_CurrentSetting_f == 1) {
      /* Sum: '<S184>/Sum' */
      mainV03_56_B.Sum_dk = (mainV03_56_B.ProportionalGain_dp +
        mainV03_56_B.Integrator_n4) + mainV03_56_B.FilterCoefficient_dn;
      mainV03_56_B.ManualSwitch_k = mainV03_56_B.Sum_dk;
    } else {
      /* Sum: '<S180>/Sum' */
      mainV03_56_B.Sum_ne = (mainV03_56_B.ProportionalGain_lh +
        mainV03_56_B.Integrator_bx) + mainV03_56_B.FilterCoefficient_im;
      mainV03_56_B.ManualSwitch_k = mainV03_56_B.Sum_ne;
    }

    /* End of ManualSwitch: '<S165>/Manual Switch' */

    /* ManualSwitch: '<S165>/Manual Switch1' incorporates:
     *  Constant: '<S165>/Lateral Position Cmd'
     */
    if (mainV03_56_P.ManualSwitch1_CurrentSetting_h == 1) {
      mainV03_56_B.ManualSwitch1_o = mainV03_56_B.ToBodyAxes3_l[1];
    } else {
      mainV03_56_B.ManualSwitch1_o = mainV03_56_P.LateralPositionCmd_Value;
    }

    /* End of ManualSwitch: '<S165>/Manual Switch1' */
    for (rtb_Sum1_eh = 0; rtb_Sum1_eh < 3; rtb_Sum1_eh++) {
      /* Product: '<S165>/To Body Axes6' */
      mainV03_56_B.ToBodyAxes6_a[rtb_Sum1_eh] = 0.0;

      /* Product: '<S165>/To Body Axes7' */
      mainV03_56_B.ToBodyAxes7_e[rtb_Sum1_eh] = 0.0;

      /* Product: '<S165>/To Body Axes5' */
      mainV03_56_B.ToBodyAxes5_b[rtb_Sum1_eh] = 0.0;

      /* Product: '<S165>/To Body Axes6' */
      mainV03_56_B.ToBodyAxes6_a[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_lt[rtb_Sum1_eh] *
        mainV03_56_B.Sensors.X_nedMeas[0];

      /* Product: '<S165>/To Body Axes7' */
      mainV03_56_B.ToBodyAxes7_e[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_lt[rtb_Sum1_eh] *
        mainV03_56_B.Sensors.VearthMeas[0];

      /* Product: '<S165>/To Body Axes5' */
      mainV03_56_B.ToBodyAxes5_b[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_lt[rtb_Sum1_eh] *
        mainV03_56_B.Sensors.VearthMeas[0];

      /* Product: '<S165>/To Body Axes6' */
      mainV03_56_B.ToBodyAxes6_a[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_lt[rtb_Sum1_eh + 3] *
        mainV03_56_B.Sensors.X_nedMeas[1];

      /* Product: '<S165>/To Body Axes7' */
      mainV03_56_B.ToBodyAxes7_e[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_lt[rtb_Sum1_eh + 3] *
        mainV03_56_B.Sensors.VearthMeas[1];

      /* Product: '<S165>/To Body Axes5' */
      mainV03_56_B.ToBodyAxes5_b[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_lt[rtb_Sum1_eh + 3] *
        mainV03_56_B.Sensors.VearthMeas[1];

      /* Product: '<S165>/To Body Axes6' */
      mainV03_56_B.ToBodyAxes6_a[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_lt[rtb_Sum1_eh + 6] *
        mainV03_56_B.Sensors.X_nedMeas[2];

      /* Product: '<S165>/To Body Axes7' */
      mainV03_56_B.ToBodyAxes7_e[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_lt[rtb_Sum1_eh + 6] *
        mainV03_56_B.Sensors.VearthMeas[2];

      /* Product: '<S165>/To Body Axes5' */
      mainV03_56_B.ToBodyAxes5_b[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_lt[rtb_Sum1_eh + 6] *
        mainV03_56_B.Sensors.VearthMeas[2];
    }

    /* Sum: '<S165>/Sum8' */
    mainV03_56_B.Sum8_h = mainV03_56_B.ManualSwitch1_o -
      mainV03_56_B.ToBodyAxes6_a[1];

    /* Gain: '<S185>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_m = mainV03_56_P.PIDLateralPosition_P *
      mainV03_56_B.Sum8_h;

    /* Integrator: '<S185>/Integrator' */
    mainV03_56_B.Integrator_iz = mainV03_56_X.Integrator_CSTATE_hz;

    /* Gain: '<S185>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_j3 = mainV03_56_P.PIDLateralPosition_D *
      mainV03_56_B.Sum8_h;

    /* Integrator: '<S185>/Filter' */
    mainV03_56_B.Filter_j = mainV03_56_X.Filter_CSTATE_l3;

    /* Sum: '<S185>/SumD' */
    mainV03_56_B.SumD_g3 = mainV03_56_B.DerivativeGain_j3 -
      mainV03_56_B.Filter_j;

    /* Gain: '<S185>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_dj = mainV03_56_P.PIDLateralPosition_N *
      mainV03_56_B.SumD_g3;

    /* Sum: '<S185>/Sum' */
    mainV03_56_B.Sum_pnh = (mainV03_56_B.ProportionalGain_m +
      mainV03_56_B.Integrator_iz) + mainV03_56_B.FilterCoefficient_dj;

    /* Saturate: '<S165>/Saturation2' */
    if (mainV03_56_B.Sum_pnh > mainV03_56_P.Saturation2_UpperSat) {
      mainV03_56_B.Saturation2_p = mainV03_56_P.Saturation2_UpperSat;
    } else if (mainV03_56_B.Sum_pnh < mainV03_56_P.Saturation2_LowerSat) {
      mainV03_56_B.Saturation2_p = mainV03_56_P.Saturation2_LowerSat;
    } else {
      mainV03_56_B.Saturation2_p = mainV03_56_B.Sum_pnh;
    }

    /* End of Saturate: '<S165>/Saturation2' */

    /* Sum: '<S165>/Sum3' */
    mainV03_56_B.Sum3_g = mainV03_56_B.Saturation2_p -
      mainV03_56_B.ToBodyAxes7_e[1];

    /* Gain: '<S183>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_m0 = mainV03_56_P.PIDController1_P_h *
      mainV03_56_B.Sum3_g;

    /* Integrator: '<S183>/Integrator' */
    mainV03_56_B.Integrator_fay = mainV03_56_X.Integrator_CSTATE_pc;

    /* Gain: '<S183>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_je = mainV03_56_P.PIDController1_D_o *
      mainV03_56_B.Sum3_g;

    /* Integrator: '<S183>/Filter' */
    mainV03_56_B.Filter_a = mainV03_56_X.Filter_CSTATE_c0;

    /* Sum: '<S183>/SumD' */
    mainV03_56_B.SumD_nt = mainV03_56_B.DerivativeGain_je -
      mainV03_56_B.Filter_a;

    /* Gain: '<S183>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_dr = mainV03_56_P.PIDController1_N_p *
      mainV03_56_B.SumD_nt;

    /* Sum: '<S183>/Sum' */
    mainV03_56_B.Sum_fr = (mainV03_56_B.ProportionalGain_m0 +
      mainV03_56_B.Integrator_fay) + mainV03_56_B.FilterCoefficient_dr;

    /* Saturate: '<S165>/Saturation1' */
    if (mainV03_56_B.Sum_fr > mainV03_56_P.Saturation1_UpperSat) {
      mainV03_56_B.Saturation1_h = mainV03_56_P.Saturation1_UpperSat;
    } else if (mainV03_56_B.Sum_fr < mainV03_56_P.Saturation1_LowerSat) {
      mainV03_56_B.Saturation1_h = mainV03_56_P.Saturation1_LowerSat;
    } else {
      mainV03_56_B.Saturation1_h = mainV03_56_B.Sum_fr;
    }

    /* End of Saturate: '<S165>/Saturation1' */

    /* Sum: '<S165>/Sum6' */
    mainV03_56_B.Sum6_p = mainV03_56_B.Saturation1_h -
      mainV03_56_B.Sensors.EulerMeas[0];

    /* Gain: '<S191>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_iq = mainV03_56_P.ProportionalRoll_P *
      mainV03_56_B.Sum6_p;

    /* Integrator: '<S191>/Integrator' */
    mainV03_56_B.Integrator_k = mainV03_56_X.Integrator_CSTATE_l2;

    /* Gain: '<S191>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_a1 = mainV03_56_P.ProportionalRoll_D *
      mainV03_56_B.Sum6_p;

    /* Integrator: '<S191>/Filter' */
    mainV03_56_B.Filter_fd = mainV03_56_X.Filter_CSTATE_a5;

    /* Sum: '<S191>/SumD' */
    mainV03_56_B.SumD_e0 = mainV03_56_B.DerivativeGain_a1 -
      mainV03_56_B.Filter_fd;

    /* Gain: '<S191>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_h = mainV03_56_P.ProportionalRoll_N *
      mainV03_56_B.SumD_e0;

    /* Sum: '<S191>/Sum' */
    mainV03_56_B.Sum_a2b = (mainV03_56_B.ProportionalGain_iq +
      mainV03_56_B.Integrator_k) + mainV03_56_B.FilterCoefficient_h;

    /* Sum: '<S165>/Sum7' */
    mainV03_56_B.Sum7_l = mainV03_56_B.Sum_a2b -
      mainV03_56_B.Sensors.OmegaMeas_body.pMeas;

    /* Gain: '<S187>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_dv = mainV03_56_P.PIDRollRate_P *
      mainV03_56_B.Sum7_l;

    /* Integrator: '<S187>/Integrator' */
    mainV03_56_B.Integrator_af = mainV03_56_X.Integrator_CSTATE_am;

    /* Gain: '<S187>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_na = mainV03_56_P.PIDRollRate_D *
      mainV03_56_B.Sum7_l;

    /* Integrator: '<S187>/Filter' */
    mainV03_56_B.Filter_k = mainV03_56_X.Filter_CSTATE_mg;

    /* Sum: '<S187>/SumD' */
    mainV03_56_B.SumD_m = mainV03_56_B.DerivativeGain_na - mainV03_56_B.Filter_k;

    /* Gain: '<S187>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_dv = mainV03_56_P.PIDRollRate_N *
      mainV03_56_B.SumD_m;

    /* Sum: '<S187>/Sum' */
    mainV03_56_B.Sum_af = (mainV03_56_B.ProportionalGain_dv +
      mainV03_56_B.Integrator_af) + mainV03_56_B.FilterCoefficient_dv;

    /* Gain: '<S197>/dRoll' */
    mainV03_56_B.dRoll_m = mainV03_56_P.dRoll_Gain * mainV03_56_B.Sum_af;

    /* Gain: '<S197>/Gain' */
    mainV03_56_B.Gain_h = mainV03_56_P.Gain_Gain_di * mainV03_56_B.dRoll_m;

    /* Gain: '<S195>/dThrottle' */
    mainV03_56_B.dThrottle_g = mainV03_56_P.dThrottle_Gain *
      mainV03_56_B.ManualSwitch_k;

    /* Bias: '<S195>/Bias' */
    mainV03_56_B.Bias_f = mainV03_56_B.dThrottle_g + mainV03_56_P.Bias_Bias;

    /* ManualSwitch: '<S165>/Manual Switch6' incorporates:
     *  Constant: '<S165>/Constant4'
     */
    if (mainV03_56_P.ManualSwitch6_CurrentSetting_n == 1) {
      mainV03_56_B.ManualSwitch6_e = mainV03_56_P.Constant4_Value_d;
    } else {
      mainV03_56_B.ManualSwitch6_e = mainV03_56_B.ManualSwitch.yaw_cmd;
    }

    /* End of ManualSwitch: '<S165>/Manual Switch6' */

    /* Sum: '<S165>/Sum9' */
    mainV03_56_B.Sum9_e = mainV03_56_B.ManualSwitch6_e -
      mainV03_56_B.Sensors.EulerMeas[2];

    /* Gain: '<S192>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_nl = mainV03_56_P.ProportionalYaw1_P *
      mainV03_56_B.Sum9_e;

    /* Integrator: '<S192>/Integrator' */
    mainV03_56_B.Integrator_ar = mainV03_56_X.Integrator_CSTATE_dkb;

    /* Gain: '<S192>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_oq = mainV03_56_P.ProportionalYaw1_D *
      mainV03_56_B.Sum9_e;

    /* Integrator: '<S192>/Filter' */
    mainV03_56_B.Filter_b = mainV03_56_X.Filter_CSTATE_oz;

    /* Sum: '<S192>/SumD' */
    mainV03_56_B.SumD_c2 = mainV03_56_B.DerivativeGain_oq -
      mainV03_56_B.Filter_b;

    /* Gain: '<S192>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_cz = mainV03_56_P.ProportionalYaw1_N *
      mainV03_56_B.SumD_c2;

    /* Sum: '<S192>/Sum' */
    mainV03_56_B.Sum_i4 = (mainV03_56_B.ProportionalGain_nl +
      mainV03_56_B.Integrator_ar) + mainV03_56_B.FilterCoefficient_cz;

    /* Sum: '<S165>/Sum10' */
    mainV03_56_B.Sum10_n = mainV03_56_B.Sum_i4 -
      mainV03_56_B.Sensors.OmegaMeas_body.rMeas;

    /* Gain: '<S188>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_p = mainV03_56_P.PIDYawRate1_P *
      mainV03_56_B.Sum10_n;

    /* Integrator: '<S188>/Integrator' */
    mainV03_56_B.Integrator_ly = mainV03_56_X.Integrator_CSTATE_i5;

    /* Gain: '<S188>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_li = mainV03_56_P.PIDYawRate1_D *
      mainV03_56_B.Sum10_n;

    /* Integrator: '<S188>/Filter' */
    mainV03_56_B.Filter_p0 = mainV03_56_X.Filter_CSTATE_px;

    /* Sum: '<S188>/SumD' */
    mainV03_56_B.SumD_lv = mainV03_56_B.DerivativeGain_li -
      mainV03_56_B.Filter_p0;

    /* Gain: '<S188>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_ix = mainV03_56_P.PIDYawRate1_N *
      mainV03_56_B.SumD_lv;

    /* Sum: '<S188>/Sum' */
    mainV03_56_B.Sum_dv = (mainV03_56_B.ProportionalGain_p +
      mainV03_56_B.Integrator_ly) + mainV03_56_B.FilterCoefficient_ix;

    /* Gain: '<S198>/dYaw' */
    mainV03_56_B.dYaw_a = mainV03_56_P.dYaw_Gain * mainV03_56_B.Sum_dv;

    /* Saturate: '<S165>/Saturation7' */
    if (mainV03_56_B.Sum_pnb > mainV03_56_P.Saturation7_UpperSat) {
      mainV03_56_B.Saturation7_k = mainV03_56_P.Saturation7_UpperSat;
    } else if (mainV03_56_B.Sum_pnb < mainV03_56_P.Saturation7_LowerSat) {
      mainV03_56_B.Saturation7_k = mainV03_56_P.Saturation7_LowerSat;
    } else {
      mainV03_56_B.Saturation7_k = mainV03_56_B.Sum_pnb;
    }

    /* End of Saturate: '<S165>/Saturation7' */

    /* Sum: '<S165>/Sum15' */
    mainV03_56_B.Sum15_g = mainV03_56_B.Saturation7_k -
      mainV03_56_B.ToBodyAxes5_b[0];

    /* Gain: '<S182>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_fn = mainV03_56_P.PIDController_P_n *
      mainV03_56_B.Sum15_g;

    /* Integrator: '<S182>/Integrator' */
    mainV03_56_B.Integrator_ne = mainV03_56_X.Integrator_CSTATE_ng;

    /* Gain: '<S182>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_fe = mainV03_56_P.PIDController_D_n *
      mainV03_56_B.Sum15_g;

    /* Integrator: '<S182>/Filter' */
    mainV03_56_B.Filter_ik = mainV03_56_X.Filter_CSTATE_d;

    /* Sum: '<S182>/SumD' */
    mainV03_56_B.SumD_mm = mainV03_56_B.DerivativeGain_fe -
      mainV03_56_B.Filter_ik;

    /* Gain: '<S182>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_g = mainV03_56_P.PIDController_N_d *
      mainV03_56_B.SumD_mm;

    /* Sum: '<S182>/Sum' */
    mainV03_56_B.Sum_ib = (mainV03_56_B.ProportionalGain_fn +
      mainV03_56_B.Integrator_ne) + mainV03_56_B.FilterCoefficient_g;

    /* Saturate: '<S165>/Saturation6' */
    if (mainV03_56_B.Sum_ib > mainV03_56_P.Saturation6_UpperSat) {
      mainV03_56_B.Saturation6_k = mainV03_56_P.Saturation6_UpperSat;
    } else if (mainV03_56_B.Sum_ib < mainV03_56_P.Saturation6_LowerSat) {
      mainV03_56_B.Saturation6_k = mainV03_56_P.Saturation6_LowerSat;
    } else {
      mainV03_56_B.Saturation6_k = mainV03_56_B.Sum_ib;
    }

    /* End of Saturate: '<S165>/Saturation6' */

    /* Sum: '<S165>/Sum12' */
    mainV03_56_B.Sum12_e = mainV03_56_B.Saturation6_k -
      mainV03_56_B.Sensors.EulerMeas[1];

    /* Gain: '<S190>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_dx = mainV03_56_P.ProportionalPitch2_P *
      mainV03_56_B.Sum12_e;

    /* Integrator: '<S190>/Integrator' */
    mainV03_56_B.Integrator_nen = mainV03_56_X.Integrator_CSTATE_cp;

    /* Gain: '<S190>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_no = mainV03_56_P.ProportionalPitch2_D *
      mainV03_56_B.Sum12_e;

    /* Integrator: '<S190>/Filter' */
    mainV03_56_B.Filter_h5 = mainV03_56_X.Filter_CSTATE_gi;

    /* Sum: '<S190>/SumD' */
    mainV03_56_B.SumD_fi = mainV03_56_B.DerivativeGain_no -
      mainV03_56_B.Filter_h5;

    /* Gain: '<S190>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_dc = mainV03_56_P.ProportionalPitch2_N *
      mainV03_56_B.SumD_fi;

    /* Sum: '<S190>/Sum' */
    mainV03_56_B.Sum_cy = (mainV03_56_B.ProportionalGain_dx +
      mainV03_56_B.Integrator_nen) + mainV03_56_B.FilterCoefficient_dc;

    /* Sum: '<S165>/Sum13' */
    mainV03_56_B.PitchRateIn_c = mainV03_56_B.Sum_cy -
      mainV03_56_B.Sensors.OmegaMeas_body.qMeas;

    /* Gain: '<S186>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_ab = mainV03_56_P.PIDPitchRate2_P *
      mainV03_56_B.PitchRateIn_c;

    /* Integrator: '<S186>/Integrator' */
    mainV03_56_B.Integrator_hd = mainV03_56_X.Integrator_CSTATE_e1;

    /* Gain: '<S186>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_cp = mainV03_56_P.PIDPitchRate2_D *
      mainV03_56_B.PitchRateIn_c;

    /* Integrator: '<S186>/Filter' */
    mainV03_56_B.Filter_dy = mainV03_56_X.Filter_CSTATE_pc;

    /* Sum: '<S186>/SumD' */
    mainV03_56_B.SumD_ei = mainV03_56_B.DerivativeGain_cp -
      mainV03_56_B.Filter_dy;

    /* Gain: '<S186>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_mr = mainV03_56_P.PIDPitchRate2_N *
      mainV03_56_B.SumD_ei;

    /* Sum: '<S186>/Sum' */
    mainV03_56_B.Sum_pf = (mainV03_56_B.ProportionalGain_ab +
      mainV03_56_B.Integrator_hd) + mainV03_56_B.FilterCoefficient_mr;

    /* Gain: '<S196>/dPitch' */
    mainV03_56_B.dPitch_h = mainV03_56_P.dPitch_Gain * mainV03_56_B.Sum_pf;

    /* Sum: '<S193>/Sum2' */
    mainV03_56_B.Throttle2noSaturation_l = ((mainV03_56_B.Gain_h +
      mainV03_56_B.Bias_f) + mainV03_56_B.dYaw_a) + mainV03_56_B.dPitch_h;

    /* Saturate: '<S193>/Saturation' */
    if (mainV03_56_B.Throttle2noSaturation_l >
        mainV03_56_P.Saturation_UpperSat_p) {
      mainV03_56_B.BusCreator1_h.Throttle2 = mainV03_56_P.Saturation_UpperSat_p;
    } else if (mainV03_56_B.Throttle2noSaturation_l <
               mainV03_56_P.Saturation_LowerSat_c) {
      mainV03_56_B.BusCreator1_h.Throttle2 = mainV03_56_P.Saturation_LowerSat_c;
    } else {
      mainV03_56_B.BusCreator1_h.Throttle2 =
        mainV03_56_B.Throttle2noSaturation_l;
    }

    /* End of Saturate: '<S193>/Saturation' */
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* Product: '<S196>/Divide' incorporates:
       *  Constant: '<S196>/Constant1'
       *  Constant: '<S196>/Constant2'
       *  Constant: '<S196>/motor2 position '
       *  Constant: '<S196>/motor2 position 1'
       *  Sum: '<S196>/Sum1'
       *  Sum: '<S196>/Sum2'
       */
      mainV03_56_B.Divide_h = (mainV03_56_P.motor2position_Value[0] -
        mainV03_56_P.Constant1_Value_n[0]) /
        (mainV03_56_P.motor2position1_Value[0] - mainV03_56_P.Constant2_Value_d
         [0]);
    }

    /* Product: '<S196>/Product' */
    mainV03_56_B.Product_ma = mainV03_56_B.dPitch_h * mainV03_56_B.Divide_h;

    /* Gain: '<S198>/Gain' */
    mainV03_56_B.Gain_ac = mainV03_56_P.Gain_Gain_g * mainV03_56_B.dYaw_a;

    /* Gain: '<S197>/Gain1' */
    mainV03_56_B.Gain1_h = mainV03_56_P.Gain1_Gain_e4 * mainV03_56_B.dRoll_m;

    /* Sum: '<S193>/Sum1' */
    mainV03_56_B.Throttle3noSaturation_d = ((mainV03_56_B.Product_ma +
      mainV03_56_B.Bias_f) + mainV03_56_B.Gain_ac) + mainV03_56_B.Gain1_h;

    /* Saturate: '<S193>/Saturation1' */
    if (mainV03_56_B.Throttle3noSaturation_d >
        mainV03_56_P.Saturation1_UpperSat_h) {
      mainV03_56_B.BusCreator1_h.Throttle3 = mainV03_56_P.Saturation1_UpperSat_h;
    } else if (mainV03_56_B.Throttle3noSaturation_d <
               mainV03_56_P.Saturation1_LowerSat_i) {
      mainV03_56_B.BusCreator1_h.Throttle3 = mainV03_56_P.Saturation1_LowerSat_i;
    } else {
      mainV03_56_B.BusCreator1_h.Throttle3 =
        mainV03_56_B.Throttle3noSaturation_d;
    }

    /* End of Saturate: '<S193>/Saturation1' */

    /* Product: '<S196>/Product1' */
    mainV03_56_B.Product1_l = mainV03_56_B.dPitch_h * mainV03_56_B.Divide_h;

    /* Sum: '<S193>/Sum3' */
    mainV03_56_B.Throttle4noSaturation_d = ((mainV03_56_B.dRoll_m +
      mainV03_56_B.Bias_f) + mainV03_56_B.dYaw_a) + mainV03_56_B.Product1_l;

    /* Saturate: '<S193>/Saturation2' */
    if (mainV03_56_B.Throttle4noSaturation_d >
        mainV03_56_P.Saturation2_UpperSat_m) {
      mainV03_56_B.BusCreator1_h.Throttle4 = mainV03_56_P.Saturation2_UpperSat_m;
    } else if (mainV03_56_B.Throttle4noSaturation_d <
               mainV03_56_P.Saturation2_LowerSat_k) {
      mainV03_56_B.BusCreator1_h.Throttle4 = mainV03_56_P.Saturation2_LowerSat_k;
    } else {
      mainV03_56_B.BusCreator1_h.Throttle4 =
        mainV03_56_B.Throttle4noSaturation_d;
    }

    /* End of Saturate: '<S193>/Saturation2' */

    /* Gain: '<S198>/Gain1' */
    mainV03_56_B.Gain1_ne = mainV03_56_P.Gain1_Gain_k * mainV03_56_B.dYaw_a;

    /* Sum: '<S193>/Sum' */
    mainV03_56_B.Throttle5noSaturation_g = ((mainV03_56_B.dPitch_h +
      mainV03_56_B.Bias_f) + mainV03_56_B.Gain1_ne) + mainV03_56_B.dRoll_m;

    /* Saturate: '<S193>/Saturation3' */
    if (mainV03_56_B.Throttle5noSaturation_g >
        mainV03_56_P.Saturation3_UpperSat_p) {
      mainV03_56_B.BusCreator1_h.Throttle5 = mainV03_56_P.Saturation3_UpperSat_p;
    } else if (mainV03_56_B.Throttle5noSaturation_g <
               mainV03_56_P.Saturation3_LowerSat_c) {
      mainV03_56_B.BusCreator1_h.Throttle5 = mainV03_56_P.Saturation3_LowerSat_c;
    } else {
      mainV03_56_B.BusCreator1_h.Throttle5 =
        mainV03_56_B.Throttle5noSaturation_g;
    }

    /* End of Saturate: '<S193>/Saturation3' */

    /* BusCreator: '<S193>/Bus Creator1' incorporates:
     *  Constant: '<S193>/Constant13'
     */
    mainV03_56_B.BusCreator1_h.Throttle1 = mainV03_56_P.Constant13_Value_j;
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* ManualSwitch: '<S165>/Manual Switch5' */
      if (mainV03_56_P.ManualSwitch5_CurrentSetting_d == 1) {
        mainV03_56_B.ManualSwitch5_n.Throttle1 = rtb_sincos_o1_idx_0;
        mainV03_56_B.ManualSwitch5_n.Throttle2 = rtb_sincos_o2_g_idx_0;
        mainV03_56_B.ManualSwitch5_n.Throttle3 = rtb_sincos_o1_idx_1;
        mainV03_56_B.ManualSwitch5_n.Throttle4 = rtb_sincos_o2_g_idx_1;
        mainV03_56_B.ManualSwitch5_n.Throttle5 = rtb_sincos_o1_idx_2;
      } else {
        mainV03_56_B.ManualSwitch5_n.Throttle1 = rtb_sincos_o2_g_idx_2;
        mainV03_56_B.ManualSwitch5_n.Throttle2 = rtb_BusCreator3_n_Throttle2;
        mainV03_56_B.ManualSwitch5_n.Throttle3 = rtb_BusCreator3_n_Throttle3;
        mainV03_56_B.ManualSwitch5_n.Throttle4 = rtb_BusCreator3_n_Throttle4;
        mainV03_56_B.ManualSwitch5_n.Throttle5 = rtb_BusCreator3_n_Throttle5;
      }

      /* End of ManualSwitch: '<S165>/Manual Switch5' */
    }

    /* ManualSwitch: '<S165>/Manual Switch3' */
    if (mainV03_56_P.ManualSwitch3_CurrentSetting_i == 1) {
      mainV03_56_B.ManualSwitch3_i = mainV03_56_B.BusCreator1_h;
    } else {
      mainV03_56_B.ManualSwitch3_i = mainV03_56_B.ManualSwitch5_n;
    }

    /* End of ManualSwitch: '<S165>/Manual Switch3' */

    /* Gain: '<S180>/Integral Gain' */
    mainV03_56_B.IntegralGain_nx = mainV03_56_P.PIDAltitudeAcceleration_I *
      mainV03_56_B.Sum11_j;

    /* Sum: '<S165>/Sum' */
    mainV03_56_B.Sum_ly = mainV03_56_B.ManualSwitch.altitude_cmd -
      mainV03_56_B.Sensors.LLAMeas.AltitudeMeas_m;

    /* Gain: '<S189>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_nb = mainV03_56_P.ProportionalAltitude_P *
      mainV03_56_B.Sum_ly;

    /* Integrator: '<S189>/Integrator' */
    mainV03_56_B.Integrator_gb = mainV03_56_X.Integrator_CSTATE_gw;

    /* Gain: '<S189>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_j33 = mainV03_56_P.ProportionalAltitude_D *
      mainV03_56_B.Sum_ly;

    /* Integrator: '<S189>/Filter' */
    mainV03_56_B.Filter_kh = mainV03_56_X.Filter_CSTATE_py;

    /* Sum: '<S189>/SumD' */
    mainV03_56_B.SumD_c3 = mainV03_56_B.DerivativeGain_j33 -
      mainV03_56_B.Filter_kh;

    /* Gain: '<S189>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_ko = mainV03_56_P.ProportionalAltitude_N *
      mainV03_56_B.SumD_c3;

    /* Sum: '<S189>/Sum' */
    mainV03_56_B.Sum_kb = (mainV03_56_B.ProportionalGain_nb +
      mainV03_56_B.Integrator_gb) + mainV03_56_B.FilterCoefficient_ko;

    /* Saturate: '<S189>/Saturate' */
    if (mainV03_56_B.Sum_kb >
        mainV03_56_P.ProportionalAltitude_UpperSaturationLimit) {
      mainV03_56_B.Saturate_iq =
        mainV03_56_P.ProportionalAltitude_UpperSaturationLimit;
    } else if (mainV03_56_B.Sum_kb <
               mainV03_56_P.ProportionalAltitude_LowerSaturationLimit) {
      mainV03_56_B.Saturate_iq =
        mainV03_56_P.ProportionalAltitude_LowerSaturationLimit;
    } else {
      mainV03_56_B.Saturate_iq = mainV03_56_B.Sum_kb;
    }

    /* End of Saturate: '<S189>/Saturate' */

    /* Sum: '<S165>/Sum1' incorporates:
     *  UnaryMinus: '<S165>/Unary Minus'
     */
    mainV03_56_B.Sum1_c = mainV03_56_B.Saturate_iq -
      (-mainV03_56_B.Sensors.VearthMeas[2]);

    /* Gain: '<S181>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_inf = mainV03_56_P.PIDAltitudeRate_D *
      mainV03_56_B.Sum1_c;

    /* Integrator: '<S181>/Filter' */
    mainV03_56_B.Filter_gp = mainV03_56_X.Filter_CSTATE_ib;

    /* Sum: '<S181>/SumD' */
    mainV03_56_B.SumD_l3 = mainV03_56_B.DerivativeGain_inf -
      mainV03_56_B.Filter_gp;

    /* Gain: '<S181>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_lz = mainV03_56_P.PIDAltitudeRate_N *
      mainV03_56_B.SumD_l3;

    /* Gain: '<S181>/Integral Gain' */
    mainV03_56_B.IntegralGain_ke = mainV03_56_P.PIDAltitudeRate_I *
      mainV03_56_B.Sum1_c;

    /* Gain: '<S182>/Integral Gain' */
    mainV03_56_B.IntegralGain_l = mainV03_56_P.PIDController_I_f *
      mainV03_56_B.Sum15_g;

    /* Gain: '<S183>/Integral Gain' */
    mainV03_56_B.IntegralGain_g4 = mainV03_56_P.PIDController1_I_l *
      mainV03_56_B.Sum3_g;

    /* Gain: '<S184>/Integral Gain' */
    mainV03_56_B.IntegralGain_hx = mainV03_56_P.PIDHeightController1_I *
      mainV03_56_B.Saturation9_k;

    /* Gain: '<S185>/Integral Gain' */
    mainV03_56_B.IntegralGain_c = mainV03_56_P.PIDLateralPosition_I *
      mainV03_56_B.Sum8_h;

    /* Gain: '<S186>/Integral Gain' */
    mainV03_56_B.IntegralGain_ho = mainV03_56_P.PIDPitchRate2_I *
      mainV03_56_B.PitchRateIn_c;

    /* Gain: '<S187>/Integral Gain' */
    mainV03_56_B.IntegralGain_kc = mainV03_56_P.PIDRollRate_I *
      mainV03_56_B.Sum7_l;

    /* Gain: '<S188>/Integral Gain' */
    mainV03_56_B.IntegralGain_m1 = mainV03_56_P.PIDYawRate1_I *
      mainV03_56_B.Sum10_n;

    /* Gain: '<S189>/Integral Gain' */
    mainV03_56_B.IntegralGain_df = mainV03_56_P.ProportionalAltitude_I *
      mainV03_56_B.Sum_ly;

    /* Gain: '<S190>/Integral Gain' */
    mainV03_56_B.IntegralGain_an = mainV03_56_P.ProportionalPitch2_I *
      mainV03_56_B.Sum12_e;

    /* Gain: '<S191>/Integral Gain' */
    mainV03_56_B.IntegralGain_pv = mainV03_56_P.ProportionalRoll_I *
      mainV03_56_B.Sum6_p;

    /* Gain: '<S192>/Integral Gain' */
    mainV03_56_B.IntegralGain_pn = mainV03_56_P.ProportionalYaw1_I *
      mainV03_56_B.Sum9_e;
  }

  /* End of Outputs for SubSystem: '<S3>/Quadcopter' */

  /* Product: '<S3>/Product7' */
  mainV03_56_B.Product7[0] = (real_T)mainV03_56_B.Compare *
    mainV03_56_B.actuators_h.deltae;
  mainV03_56_B.Product7[1] = (real_T)mainV03_56_B.Compare *
    mainV03_56_B.actuators_h.deltar;
  mainV03_56_B.Product7[2] = (real_T)mainV03_56_B.Compare *
    mainV03_56_B.actuators_h.deltafr;
  mainV03_56_B.Product7[3] = (real_T)mainV03_56_B.Compare *
    mainV03_56_B.actuators_h.deltafl;

  /* RelationalOperator: '<S229>/Compare' incorporates:
   *  Constant: '<S229>/Constant'
   */
  mainV03_56_B.Compare_j = (mainV03_56_B.DataTypeConversion ==
    mainV03_56_P.CompareToConstant1_const);
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* SignalConversion: '<S3>/HiddenBuf_InsertedFor_Quadcopter --> Fixed-Wing_at_inport_2' */
    mainV03_56_B.HiddenBuf_InsertedFor_QuadcopterFixedWing_at_inport_2 =
      mainV03_56_B.Compare_j;

    /* Outputs for Enabled SubSystem: '<S3>/Quadcopter --> Fixed-Wing' incorporates:
     *  EnablePort: '<S166>/Enable'
     */
    if (rtmIsMajorTimeStep(mainV03_56_M)) {
      if (mainV03_56_B.HiddenBuf_InsertedFor_QuadcopterFixedWing_at_inport_2) {
        if (!mainV03_56_DW.QuadcopterFixedWing_MODE) {
          mainV03_56_DW.QuadcopterFixedWing_MODE = true;
        }
      } else {
        if (mainV03_56_DW.QuadcopterFixedWing_MODE) {
          mainV03_56_DW.QuadcopterFixedWing_MODE = false;
        }
      }
    }

    /* End of Outputs for SubSystem: '<S3>/Quadcopter --> Fixed-Wing' */
  }

  /* Outputs for Enabled SubSystem: '<S3>/Quadcopter --> Fixed-Wing' incorporates:
   *  EnablePort: '<S166>/Enable'
   */
  if (mainV03_56_DW.QuadcopterFixedWing_MODE) {
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* BusCreator: '<S166>/Bus Creator1' incorporates:
       *  Constant: '<S166>/Constant20'
       *  Constant: '<S166>/Constant21'
       *  Constant: '<S166>/Constant22'
       *  Constant: '<S166>/Constant23'
       *  Constant: '<S166>/Constant24'
       */
      rtb_BusCreator1_Throttle1 = mainV03_56_P.Constant21_Value_fz;
      rtb_BusCreator1_Throttle2 = mainV03_56_P.Constant22_Value_f;
      rtb_BusCreator1_Throttle3 = mainV03_56_P.Constant23_Value_f;
      rtb_BusCreator1_Throttle4 = mainV03_56_P.Constant24_Value_e;
      rtb_BusCreator1_Throttle5 = mainV03_56_P.Constant20_Value_e;
    }

    /* SignalConversion: '<S222>/TmpSignal ConversionAtsincosInport1' incorporates:
     *  Constant: '<S166>/Constant6'
     */
    mainV03_56_B.TmpSignalConversionAtsincosInport1_f[0] =
      mainV03_56_B.Sensors.EulerMeas[2];
    mainV03_56_B.TmpSignalConversionAtsincosInport1_f[1] =
      mainV03_56_P.Constant6_Value_j[0];
    mainV03_56_B.TmpSignalConversionAtsincosInport1_f[2] =
      mainV03_56_P.Constant6_Value_j[1];

    /* Trigonometry: '<S222>/sincos' */
    rtb_sincos_o1_idx_0 = sin(mainV03_56_B.TmpSignalConversionAtsincosInport1_f
      [0]);
    rtb_sincos_o2_g_idx_0 = cos
      (mainV03_56_B.TmpSignalConversionAtsincosInport1_f[0]);
    rtb_sincos_o1_idx_1 = sin(mainV03_56_B.TmpSignalConversionAtsincosInport1_f
      [1]);
    rtb_sincos_o2_g_idx_1 = cos
      (mainV03_56_B.TmpSignalConversionAtsincosInport1_f[1]);
    rtb_sincos_o1_idx_2 = sin(mainV03_56_B.TmpSignalConversionAtsincosInport1_f
      [2]);
    rtb_sincos_o2_g_idx_2 = cos
      (mainV03_56_B.TmpSignalConversionAtsincosInport1_f[2]);

    /* Fcn: '<S222>/Fcn11' */
    mainV03_56_B.VectorConcatenate_e[0] = rtb_sincos_o2_g_idx_1 *
      rtb_sincos_o2_g_idx_0;

    /* Fcn: '<S222>/Fcn21' */
    mainV03_56_B.VectorConcatenate_e[1] = rtb_sincos_o1_idx_2 *
      rtb_sincos_o1_idx_1 * rtb_sincos_o2_g_idx_0 - rtb_sincos_o2_g_idx_2 *
      rtb_sincos_o1_idx_0;

    /* Fcn: '<S222>/Fcn31' */
    mainV03_56_B.VectorConcatenate_e[2] = rtb_sincos_o2_g_idx_2 *
      rtb_sincos_o1_idx_1 * rtb_sincos_o2_g_idx_0 + rtb_sincos_o1_idx_2 *
      rtb_sincos_o1_idx_0;

    /* Fcn: '<S222>/Fcn12' */
    mainV03_56_B.VectorConcatenate_e[3] = rtb_sincos_o2_g_idx_1 *
      rtb_sincos_o1_idx_0;

    /* Fcn: '<S222>/Fcn22' */
    mainV03_56_B.VectorConcatenate_e[4] = rtb_sincos_o1_idx_2 *
      rtb_sincos_o1_idx_1 * rtb_sincos_o1_idx_0 + rtb_sincos_o2_g_idx_2 *
      rtb_sincos_o2_g_idx_0;

    /* Fcn: '<S222>/Fcn32' */
    mainV03_56_B.VectorConcatenate_e[5] = rtb_sincos_o2_g_idx_2 *
      rtb_sincos_o1_idx_1 * rtb_sincos_o1_idx_0 - rtb_sincos_o1_idx_2 *
      rtb_sincos_o2_g_idx_0;

    /* Fcn: '<S222>/Fcn13' */
    mainV03_56_B.VectorConcatenate_e[6] = -rtb_sincos_o1_idx_1;

    /* Fcn: '<S222>/Fcn23' */
    mainV03_56_B.VectorConcatenate_e[7] = rtb_sincos_o1_idx_2 *
      rtb_sincos_o2_g_idx_1;

    /* Fcn: '<S222>/Fcn33' */
    mainV03_56_B.VectorConcatenate_e[8] = rtb_sincos_o2_g_idx_2 *
      rtb_sincos_o2_g_idx_1;

    /* Product: '<S166>/To Body Axes3' */
    for (rtb_Sum1_eh = 0; rtb_Sum1_eh < 3; rtb_Sum1_eh++) {
      mainV03_56_B.ToBodyAxes3[rtb_Sum1_eh] = 0.0;
      mainV03_56_B.ToBodyAxes3[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_e[rtb_Sum1_eh] *
        mainV03_56_B.ManualSwitch.Position_cmd[0];
      mainV03_56_B.ToBodyAxes3[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_e[rtb_Sum1_eh + 3] *
        mainV03_56_B.ManualSwitch.Position_cmd[1];
      mainV03_56_B.ToBodyAxes3[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_e[rtb_Sum1_eh + 6] *
        mainV03_56_B.ManualSwitch.Position_cmd[2];
    }

    /* End of Product: '<S166>/To Body Axes3' */

    /* ManualSwitch: '<S166>/Manual Switch10' incorporates:
     *  Constant: '<S166>/Constant18'
     */
    if (mainV03_56_P.ManualSwitch10_CurrentSetting == 1) {
      mainV03_56_B.ManualSwitch10 = mainV03_56_P.Constant18_Value_e;
    } else {
      mainV03_56_B.ManualSwitch10 = mainV03_56_B.ToBodyAxes3[2];
    }

    /* End of ManualSwitch: '<S166>/Manual Switch10' */

    /* Sum: '<S166>/Sum18' */
    mainV03_56_B.Sum18 = mainV03_56_B.ManualSwitch10 -
      mainV03_56_B.Sensors.LLAMeas.AltitudeMeas_m;

    /* Saturate: '<S166>/Saturation11' */
    if (mainV03_56_B.Sum18 > mainV03_56_P.Saturation11_UpperSat) {
      mainV03_56_B.Saturation11 = mainV03_56_P.Saturation11_UpperSat;
    } else if (mainV03_56_B.Sum18 < mainV03_56_P.Saturation11_LowerSat) {
      mainV03_56_B.Saturation11 = mainV03_56_P.Saturation11_LowerSat;
    } else {
      mainV03_56_B.Saturation11 = mainV03_56_B.Sum18;
    }

    /* End of Saturate: '<S166>/Saturation11' */

    /* Gain: '<S201>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_k = mainV03_56_P.PHeight1_P *
      mainV03_56_B.Saturation11;

    /* Integrator: '<S201>/Integrator' */
    mainV03_56_B.Integrator_gn = mainV03_56_X.Integrator_CSTATE_n1;

    /* Gain: '<S201>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_da = mainV03_56_P.PHeight1_D *
      mainV03_56_B.Saturation11;

    /* Integrator: '<S201>/Filter' */
    mainV03_56_B.Filter_fi = mainV03_56_X.Filter_CSTATE_nt;

    /* Sum: '<S201>/SumD' */
    mainV03_56_B.SumD_a = mainV03_56_B.DerivativeGain_da -
      mainV03_56_B.Filter_fi;

    /* Gain: '<S201>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_j = mainV03_56_P.PHeight1_N *
      mainV03_56_B.SumD_a;

    /* Sum: '<S201>/Sum' */
    mainV03_56_B.Sum_pru = (mainV03_56_B.ProportionalGain_k +
      mainV03_56_B.Integrator_gn) + mainV03_56_B.FilterCoefficient_j;

    /* Saturate: '<S201>/Saturate' */
    if (mainV03_56_B.Sum_pru > mainV03_56_P.PHeight1_UpperSaturationLimit) {
      mainV03_56_B.Saturate_m = mainV03_56_P.PHeight1_UpperSaturationLimit;
    } else if (mainV03_56_B.Sum_pru < mainV03_56_P.PHeight1_LowerSaturationLimit)
    {
      mainV03_56_B.Saturate_m = mainV03_56_P.PHeight1_LowerSaturationLimit;
    } else {
      mainV03_56_B.Saturate_m = mainV03_56_B.Sum_pru;
    }

    /* End of Saturate: '<S201>/Saturate' */

    /* Sum: '<S166>/Sum19' */
    mainV03_56_B.Sum19 = mainV03_56_B.Saturate_m -
      mainV03_56_B.Sensors.alphaMeas;

    /* Gain: '<S207>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_kl = mainV03_56_P.PIDController5_P_n *
      mainV03_56_B.Sum19;

    /* Integrator: '<S207>/Integrator' */
    mainV03_56_B.Integrator_i = mainV03_56_X.Integrator_CSTATE_e;

    /* Gain: '<S207>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_g = mainV03_56_P.PIDController5_D_f *
      mainV03_56_B.Sum19;

    /* Integrator: '<S207>/Filter' */
    mainV03_56_B.Filter_o = mainV03_56_X.Filter_CSTATE_m;

    /* Sum: '<S207>/SumD' */
    mainV03_56_B.SumD_ie = mainV03_56_B.DerivativeGain_g - mainV03_56_B.Filter_o;

    /* Gain: '<S207>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_f = mainV03_56_P.PIDController5_N_p *
      mainV03_56_B.SumD_ie;

    /* Sum: '<S207>/Sum' */
    mainV03_56_B.Sum_a2 = (mainV03_56_B.ProportionalGain_kl +
      mainV03_56_B.Integrator_i) + mainV03_56_B.FilterCoefficient_f;

    /* Sum: '<S166>/Sum20' */
    mainV03_56_B.PitchRateIn_j = mainV03_56_B.Sum_a2 -
      mainV03_56_B.Sensors.OmegaMeas_body.qMeas;

    /* Gain: '<S206>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_nu = mainV03_56_P.PIDController4_P_k *
      mainV03_56_B.PitchRateIn_j;

    /* Integrator: '<S206>/Integrator' */
    mainV03_56_B.Integrator_n = mainV03_56_X.Integrator_CSTATE_il;

    /* Gain: '<S206>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_i = mainV03_56_P.PIDController4_D_n *
      mainV03_56_B.PitchRateIn_j;

    /* Integrator: '<S206>/Filter' */
    mainV03_56_B.Filter_e = mainV03_56_X.Filter_CSTATE_i;

    /* Sum: '<S206>/SumD' */
    mainV03_56_B.SumD_g = mainV03_56_B.DerivativeGain_i - mainV03_56_B.Filter_e;

    /* Gain: '<S206>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_ln = mainV03_56_P.PIDController4_N_p *
      mainV03_56_B.SumD_g;

    /* Sum: '<S206>/Sum' */
    mainV03_56_B.Sum_ph = (mainV03_56_B.ProportionalGain_nu +
      mainV03_56_B.Integrator_n) + mainV03_56_B.FilterCoefficient_ln;

    /* Saturate: '<S206>/Saturate' */
    if (mainV03_56_B.Sum_ph > mainV03_56_P.PIDController4_UpperSaturationLimit)
    {
      mainV03_56_B.Saturate_k = mainV03_56_P.PIDController4_UpperSaturationLimit;
    } else if (mainV03_56_B.Sum_ph <
               mainV03_56_P.PIDController4_LowerSaturationLimit) {
      mainV03_56_B.Saturate_k = mainV03_56_P.PIDController4_LowerSaturationLimit;
    } else {
      mainV03_56_B.Saturate_k = mainV03_56_B.Sum_ph;
    }

    /* End of Saturate: '<S206>/Saturate' */
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* Gain: '<S166>/Gain4' incorporates:
       *  Constant: '<S166>/Constant7'
       */
      mainV03_56_B.deltae_p = mainV03_56_P.Gain4_Gain_f * mainV03_56_P.delta_e;
    }

    /* ManualSwitch: '<S166>/Manual Switch2' */
    if (mainV03_56_P.ManualSwitch2_CurrentSetting_c == 1) {
      mainV03_56_B.deltae_l = mainV03_56_B.Saturate_k;
    } else {
      mainV03_56_B.deltae_l = mainV03_56_B.deltae_p;
    }

    /* End of ManualSwitch: '<S166>/Manual Switch2' */

    /* Sum: '<S166>/Sum22' incorporates:
     *  Constant: '<S166>/Constant17'
     */
    mainV03_56_B.Sum22 = mainV03_56_P.Constant17_Value_n -
      mainV03_56_B.Sensors.EulerMeas[2];

    /* Gain: '<S209>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_gi = mainV03_56_P.PIDController7_P_a *
      mainV03_56_B.Sum22;

    /* Integrator: '<S209>/Integrator' */
    mainV03_56_B.Integrator_en = mainV03_56_X.Integrator_CSTATE_mj;

    /* Gain: '<S209>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_eh = mainV03_56_P.PIDController7_D_i *
      mainV03_56_B.Sum22;

    /* Integrator: '<S209>/Filter' */
    mainV03_56_B.Filter_l = mainV03_56_X.Filter_CSTATE_kd;

    /* Sum: '<S209>/SumD' */
    mainV03_56_B.SumD_ch = mainV03_56_B.DerivativeGain_eh -
      mainV03_56_B.Filter_l;

    /* Gain: '<S209>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_mb = mainV03_56_P.PIDController7_N_j *
      mainV03_56_B.SumD_ch;

    /* Sum: '<S209>/Sum' */
    mainV03_56_B.Sum_mw = (mainV03_56_B.ProportionalGain_gi +
      mainV03_56_B.Integrator_en) + mainV03_56_B.FilterCoefficient_mb;

    /* Sum: '<S166>/Sum21' */
    mainV03_56_B.Sum21 = mainV03_56_B.Sum_mw -
      mainV03_56_B.Sensors.OmegaMeas_body.rMeas;

    /* Gain: '<S208>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_c = mainV03_56_P.PIDController6_P *
      mainV03_56_B.Sum21;

    /* Integrator: '<S208>/Integrator' */
    mainV03_56_B.Integrator_d = mainV03_56_X.Integrator_CSTATE_dk;

    /* Gain: '<S208>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_n = mainV03_56_P.PIDController6_D *
      mainV03_56_B.Sum21;

    /* Integrator: '<S208>/Filter' */
    mainV03_56_B.Filter_d2 = mainV03_56_X.Filter_CSTATE_ii;

    /* Sum: '<S208>/SumD' */
    mainV03_56_B.SumD_p = mainV03_56_B.DerivativeGain_n - mainV03_56_B.Filter_d2;

    /* Gain: '<S208>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_n = mainV03_56_P.PIDController6_N *
      mainV03_56_B.SumD_p;
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* Gain: '<S166>/Gain5' incorporates:
       *  Constant: '<S166>/Constant8'
       */
      mainV03_56_B.deltar_p = mainV03_56_P.Gain5_Gain_p *
        mainV03_56_P.Constant8_Value_l;
    }

    /* ManualSwitch: '<S166>/Manual Switch4' */
    if (mainV03_56_P.ManualSwitch4_CurrentSetting_b == 1) {
      /* Sum: '<S208>/Sum' */
      mainV03_56_B.Sum_ae = (mainV03_56_B.ProportionalGain_c +
        mainV03_56_B.Integrator_d) + mainV03_56_B.FilterCoefficient_n;
      mainV03_56_B.deltar_j = mainV03_56_B.Sum_ae;
    } else {
      mainV03_56_B.deltar_j = mainV03_56_B.deltar_p;
    }

    /* End of ManualSwitch: '<S166>/Manual Switch4' */

    /* Sum: '<S166>/Sum24' incorporates:
     *  Constant: '<S166>/Constant19'
     */
    mainV03_56_B.Sum24 = mainV03_56_P.Constant19_Value -
      mainV03_56_B.Sensors.EulerMeas[0];

    /* Gain: '<S211>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_b = mainV03_56_P.PIDController9_P *
      mainV03_56_B.Sum24;

    /* Integrator: '<S211>/Integrator' */
    mainV03_56_B.Integrator_b = mainV03_56_X.Integrator_CSTATE_l;

    /* Gain: '<S211>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_ha = mainV03_56_P.PIDController9_D *
      mainV03_56_B.Sum24;

    /* Integrator: '<S211>/Filter' */
    mainV03_56_B.Filter_i = mainV03_56_X.Filter_CSTATE_c;

    /* Sum: '<S211>/SumD' */
    mainV03_56_B.SumD_gv = mainV03_56_B.DerivativeGain_ha -
      mainV03_56_B.Filter_i;

    /* Gain: '<S211>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_nn = mainV03_56_P.PIDController9_N *
      mainV03_56_B.SumD_gv;

    /* Sum: '<S211>/Sum' */
    mainV03_56_B.Sum_bm = (mainV03_56_B.ProportionalGain_b +
      mainV03_56_B.Integrator_b) + mainV03_56_B.FilterCoefficient_nn;

    /* Sum: '<S166>/Sum23' */
    mainV03_56_B.Sum23 = mainV03_56_B.Sum_bm -
      mainV03_56_B.Sensors.OmegaMeas_body.pMeas;

    /* Gain: '<S210>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_cc = mainV03_56_P.PIDController8_P *
      mainV03_56_B.Sum23;

    /* Integrator: '<S210>/Integrator' */
    mainV03_56_B.Integrator_fy = mainV03_56_X.Integrator_CSTATE_p0;

    /* Gain: '<S210>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_c = mainV03_56_P.PIDController8_D *
      mainV03_56_B.Sum23;

    /* Integrator: '<S210>/Filter' */
    mainV03_56_B.Filter_g = mainV03_56_X.Filter_CSTATE_a;

    /* Sum: '<S210>/SumD' */
    mainV03_56_B.SumD_b = mainV03_56_B.DerivativeGain_c - mainV03_56_B.Filter_g;

    /* Gain: '<S210>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_d = mainV03_56_P.PIDController8_N *
      mainV03_56_B.SumD_b;

    /* Sum: '<S210>/Sum' */
    mainV03_56_B.Sum_gh = (mainV03_56_B.ProportionalGain_cc +
      mainV03_56_B.Integrator_fy) + mainV03_56_B.FilterCoefficient_d;
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* ManualSwitch: '<S166>/Manual Switch13' incorporates:
       *  Constant: '<S166>/Constant2'
       *  Constant: '<S166>/Constant25'
       */
      if (mainV03_56_P.ManualSwitch13_CurrentSetting == 1) {
        mainV03_56_B.ManualSwitch13 = mainV03_56_P.Constant25_Value;
      } else {
        mainV03_56_B.ManualSwitch13 = mainV03_56_P.Constant2_Value_o;
      }

      /* End of ManualSwitch: '<S166>/Manual Switch13' */

      /* Gain: '<S166>/Gain6' incorporates:
       *  Constant: '<S166>/Constant14'
       */
      mainV03_56_B.deltafr_e = mainV03_56_P.Gain6_Gain_h *
        mainV03_56_P.Constant14_Value_p;

      /* Gain: '<S166>/Gain7' incorporates:
       *  Constant: '<S166>/Constant16'
       */
      mainV03_56_B.deltafl_e = mainV03_56_P.Gain7_Gain_f *
        mainV03_56_P.Constant16_Value_n;
    }

    /* ManualSwitch: '<S166>/Manual Switch7' */
    if (mainV03_56_P.ManualSwitch7_CurrentSetting_f == 1) {
      /* Sum: '<S166>/Sum25' */
      mainV03_56_B.Sum25 = mainV03_56_B.Sum_gh + mainV03_56_B.ManualSwitch13;
      mainV03_56_B.deltafr_n = mainV03_56_B.Sum25;
    } else {
      mainV03_56_B.deltafr_n = mainV03_56_B.deltafr_e;
    }

    /* End of ManualSwitch: '<S166>/Manual Switch7' */

    /* ManualSwitch: '<S166>/Manual Switch8' */
    if (mainV03_56_P.ManualSwitch8_CurrentSetting == 1) {
      /* Gain: '<S166>/Gain8' */
      mainV03_56_B.Gain8_l = mainV03_56_P.Gain8_Gain_d * mainV03_56_B.Sum_gh;

      /* Sum: '<S166>/Sum17' */
      mainV03_56_B.Sum17 = mainV03_56_B.Gain8_l + mainV03_56_B.ManualSwitch13;
      mainV03_56_B.deltafl_f = mainV03_56_B.Sum17;
    } else {
      mainV03_56_B.deltafl_f = mainV03_56_B.deltafl_e;
    }

    /* End of ManualSwitch: '<S166>/Manual Switch8' */

    /* BusCreator: '<S166>/Bus Creator2' */
    mainV03_56_B.BusCreator2.deltae = mainV03_56_B.deltae_l;
    mainV03_56_B.BusCreator2.deltar = mainV03_56_B.deltar_j;
    mainV03_56_B.BusCreator2.deltafr = mainV03_56_B.deltafr_n;
    mainV03_56_B.BusCreator2.deltafl = mainV03_56_B.deltafl_f;
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* BusCreator: '<S166>/Bus Creator3' incorporates:
       *  Constant: '<S166>/deltat_10'
       *  Constant: '<S166>/deltat_6'
       *  Constant: '<S166>/deltat_7'
       *  Constant: '<S166>/deltat_8'
       *  Constant: '<S166>/deltat_9'
       */
      rtb_BusCreator3_o_Throttle1 = mainV03_56_P.deltat_1;
      rtb_BusCreator3_o_Throttle2 = mainV03_56_P.deltat_2;
      rtb_BusCreator3_o_Throttle3 = mainV03_56_P.deltat_3;
      rtb_BusCreator3_o_Throttle4 = mainV03_56_P.deltat_4;
      rtb_BusCreator3_o_Throttle5 = mainV03_56_P.deltat_5;
    }

    /* ManualSwitch: '<S166>/Manual Switch12' incorporates:
     *  Constant: '<S166>/Constant13'
     */
    if (mainV03_56_P.ManualSwitch12_CurrentSetting_k == 1) {
      mainV03_56_B.ManualSwitch12 = mainV03_56_B.ToBodyAxes3[0];
    } else {
      mainV03_56_B.ManualSwitch12 = mainV03_56_P.Constant13_Value_m;
    }

    /* End of ManualSwitch: '<S166>/Manual Switch12' */

    /* Product: '<S166>/To Body Axes4' */
    for (rtb_Sum1_eh = 0; rtb_Sum1_eh < 3; rtb_Sum1_eh++) {
      mainV03_56_B.ToBodyAxes4[rtb_Sum1_eh] = 0.0;
      mainV03_56_B.ToBodyAxes4[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_e[rtb_Sum1_eh] *
        mainV03_56_B.Sensors.X_nedMeas[0];
      mainV03_56_B.ToBodyAxes4[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_e[rtb_Sum1_eh + 3] *
        mainV03_56_B.Sensors.X_nedMeas[1];
      mainV03_56_B.ToBodyAxes4[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_e[rtb_Sum1_eh + 6] *
        mainV03_56_B.Sensors.X_nedMeas[2];
    }

    /* End of Product: '<S166>/To Body Axes4' */

    /* Sum: '<S166>/Sum14' */
    mainV03_56_B.Sum14 = mainV03_56_B.ManualSwitch12 - mainV03_56_B.ToBodyAxes4
      [0];

    /* Gain: '<S200>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_eb = mainV03_56_P.ESTE2_D_i * mainV03_56_B.Sum14;

    /* Integrator: '<S200>/Filter' */
    mainV03_56_B.Filter_o0 = mainV03_56_X.Filter_CSTATE_mr;

    /* Sum: '<S200>/SumD' */
    mainV03_56_B.SumD_o = mainV03_56_B.DerivativeGain_eb -
      mainV03_56_B.Filter_o0;

    /* Gain: '<S200>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_mw = mainV03_56_P.ESTE2_N_d *
      mainV03_56_B.SumD_o;

    /* Gain: '<S200>/Integral Gain' */
    mainV03_56_B.IntegralGain_h = mainV03_56_P.ESTE2_I_p * mainV03_56_B.Sum14;

    /* Integrator: '<S200>/Integrator' */
    mainV03_56_B.Integrator_m = mainV03_56_X.Integrator_CSTATE_mx;

    /* Gain: '<S200>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_kp = mainV03_56_P.ESTE2_P_p *
      mainV03_56_B.Sum14;

    /* Sum: '<S200>/Sum' */
    mainV03_56_B.Sum_oz = (mainV03_56_B.ProportionalGain_kp +
      mainV03_56_B.Integrator_m) + mainV03_56_B.FilterCoefficient_mw;

    /* Sum: '<S166>/Sum2' */
    mainV03_56_B.delta = mainV03_56_B.ToBodyAxes3[2] -
      mainV03_56_B.Sensors.LLAMeas.AltitudeMeas_m;

    /* Saturate: '<S166>/Saturation9' */
    if (mainV03_56_B.delta > mainV03_56_P.Saturation9_UpperSat_i) {
      mainV03_56_B.Saturation9 = mainV03_56_P.Saturation9_UpperSat_i;
    } else if (mainV03_56_B.delta < mainV03_56_P.Saturation9_LowerSat_j) {
      mainV03_56_B.Saturation9 = mainV03_56_P.Saturation9_LowerSat_j;
    } else {
      mainV03_56_B.Saturation9 = mainV03_56_B.delta;
    }

    /* End of Saturate: '<S166>/Saturation9' */

    /* Gain: '<S212>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_n0 = mainV03_56_P.PIDHeightController1_P_h *
      mainV03_56_B.Saturation9;

    /* Integrator: '<S212>/Integrator' */
    mainV03_56_B.Integrator_h = mainV03_56_X.Integrator_CSTATE_g;

    /* Gain: '<S212>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_b = mainV03_56_P.PIDHeightController1_D_l *
      mainV03_56_B.Saturation9;

    /* Integrator: '<S212>/Filter' */
    mainV03_56_B.Filter_n = mainV03_56_X.Filter_CSTATE_o;

    /* Sum: '<S212>/SumD' */
    mainV03_56_B.SumD_k = mainV03_56_B.DerivativeGain_b - mainV03_56_B.Filter_n;

    /* Gain: '<S212>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_p = mainV03_56_P.PIDHeightController1_N_l *
      mainV03_56_B.SumD_k;
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* Saturate: '<S166>/Saturation3' incorporates:
       *  Constant: '<S166>/Constant'
       */
      if (mainV03_56_P.Constant_Value_lj > mainV03_56_P.Saturation3_UpperSat_n)
      {
        mainV03_56_B.Saturation3 = mainV03_56_P.Saturation3_UpperSat_n;
      } else if (mainV03_56_P.Constant_Value_lj <
                 mainV03_56_P.Saturation3_LowerSat_e) {
        mainV03_56_B.Saturation3 = mainV03_56_P.Saturation3_LowerSat_e;
      } else {
        mainV03_56_B.Saturation3 = mainV03_56_P.Constant_Value_lj;
      }

      /* End of Saturate: '<S166>/Saturation3' */
    }

    /* Sum: '<S166>/Sum11' incorporates:
     *  UnaryMinus: '<S166>/Unary Minus1'
     */
    mainV03_56_B.Sum11 = mainV03_56_B.Saturation3 -
      (-mainV03_56_B.Sensors.AccelMeas_body[2]);

    /* Gain: '<S202>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_i = mainV03_56_P.PIDAltitudeAcceleration_P_f *
      mainV03_56_B.Sum11;

    /* Integrator: '<S202>/Integrator' */
    mainV03_56_B.Integrator_fa = mainV03_56_X.Integrator_CSTATE_o;

    /* Gain: '<S202>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_gz = mainV03_56_P.PIDAltitudeAcceleration_D_l *
      mainV03_56_B.Sum11;

    /* Integrator: '<S202>/Filter' */
    mainV03_56_B.Filter_gk = mainV03_56_X.Filter_CSTATE_mk;

    /* Sum: '<S202>/SumD' */
    mainV03_56_B.SumD_kz = mainV03_56_B.DerivativeGain_gz -
      mainV03_56_B.Filter_gk;

    /* Gain: '<S202>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_b = mainV03_56_P.PIDAltitudeAcceleration_N_l *
      mainV03_56_B.SumD_kz;

    /* ManualSwitch: '<S166>/Manual Switch' */
    if (mainV03_56_P.ManualSwitch_CurrentSetting_ff == 1) {
      /* Sum: '<S212>/Sum' */
      mainV03_56_B.Sum_mb = (mainV03_56_B.ProportionalGain_n0 +
        mainV03_56_B.Integrator_h) + mainV03_56_B.FilterCoefficient_p;
      mainV03_56_B.ManualSwitch_o = mainV03_56_B.Sum_mb;
    } else {
      /* Sum: '<S202>/Sum' */
      mainV03_56_B.Sum_hg = (mainV03_56_B.ProportionalGain_i +
        mainV03_56_B.Integrator_fa) + mainV03_56_B.FilterCoefficient_b;
      mainV03_56_B.ManualSwitch_o = mainV03_56_B.Sum_hg;
    }

    /* End of ManualSwitch: '<S166>/Manual Switch' */

    /* ManualSwitch: '<S166>/Manual Switch1' incorporates:
     *  Constant: '<S166>/Lateral Position Cmd'
     */
    if (mainV03_56_P.ManualSwitch1_CurrentSetting_n == 1) {
      mainV03_56_B.ManualSwitch1_i = mainV03_56_B.ToBodyAxes3[1];
    } else {
      mainV03_56_B.ManualSwitch1_i = mainV03_56_P.LateralPositionCmd_Value_c;
    }

    /* End of ManualSwitch: '<S166>/Manual Switch1' */
    for (rtb_Sum1_eh = 0; rtb_Sum1_eh < 3; rtb_Sum1_eh++) {
      /* Product: '<S166>/To Body Axes6' */
      mainV03_56_B.ToBodyAxes6[rtb_Sum1_eh] = 0.0;

      /* Product: '<S166>/To Body Axes7' */
      mainV03_56_B.ToBodyAxes7[rtb_Sum1_eh] = 0.0;

      /* Product: '<S166>/To Body Axes5' */
      mainV03_56_B.ToBodyAxes5[rtb_Sum1_eh] = 0.0;

      /* Product: '<S166>/To Body Axes6' */
      mainV03_56_B.ToBodyAxes6[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_e[rtb_Sum1_eh] *
        mainV03_56_B.Sensors.X_nedMeas[0];

      /* Product: '<S166>/To Body Axes7' */
      mainV03_56_B.ToBodyAxes7[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_e[rtb_Sum1_eh] *
        mainV03_56_B.Sensors.VearthMeas[0];

      /* Product: '<S166>/To Body Axes5' */
      mainV03_56_B.ToBodyAxes5[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_e[rtb_Sum1_eh] *
        mainV03_56_B.Sensors.VearthMeas[0];

      /* Product: '<S166>/To Body Axes6' */
      mainV03_56_B.ToBodyAxes6[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_e[rtb_Sum1_eh + 3] *
        mainV03_56_B.Sensors.X_nedMeas[1];

      /* Product: '<S166>/To Body Axes7' */
      mainV03_56_B.ToBodyAxes7[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_e[rtb_Sum1_eh + 3] *
        mainV03_56_B.Sensors.VearthMeas[1];

      /* Product: '<S166>/To Body Axes5' */
      mainV03_56_B.ToBodyAxes5[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_e[rtb_Sum1_eh + 3] *
        mainV03_56_B.Sensors.VearthMeas[1];

      /* Product: '<S166>/To Body Axes6' */
      mainV03_56_B.ToBodyAxes6[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_e[rtb_Sum1_eh + 6] *
        mainV03_56_B.Sensors.X_nedMeas[2];

      /* Product: '<S166>/To Body Axes7' */
      mainV03_56_B.ToBodyAxes7[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_e[rtb_Sum1_eh + 6] *
        mainV03_56_B.Sensors.VearthMeas[2];

      /* Product: '<S166>/To Body Axes5' */
      mainV03_56_B.ToBodyAxes5[rtb_Sum1_eh] +=
        mainV03_56_B.VectorConcatenate_e[rtb_Sum1_eh + 6] *
        mainV03_56_B.Sensors.VearthMeas[2];
    }

    /* Sum: '<S166>/Sum8' */
    mainV03_56_B.Sum8 = mainV03_56_B.ManualSwitch1_i - mainV03_56_B.ToBodyAxes6
      [1];

    /* Gain: '<S213>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_aq = mainV03_56_P.PIDLateralPosition_P_g *
      mainV03_56_B.Sum8;

    /* Integrator: '<S213>/Integrator' */
    mainV03_56_B.Integrator_j = mainV03_56_X.Integrator_CSTATE_io;

    /* Gain: '<S213>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_df = mainV03_56_P.PIDLateralPosition_D_i *
      mainV03_56_B.Sum8;

    /* Integrator: '<S213>/Filter' */
    mainV03_56_B.Filter_if = mainV03_56_X.Filter_CSTATE_p;

    /* Sum: '<S213>/SumD' */
    mainV03_56_B.SumD_ar = mainV03_56_B.DerivativeGain_df -
      mainV03_56_B.Filter_if;

    /* Gain: '<S213>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_kl = mainV03_56_P.PIDLateralPosition_N_p *
      mainV03_56_B.SumD_ar;

    /* Sum: '<S213>/Sum' */
    mainV03_56_B.Sum_mf = (mainV03_56_B.ProportionalGain_aq +
      mainV03_56_B.Integrator_j) + mainV03_56_B.FilterCoefficient_kl;

    /* Saturate: '<S166>/Saturation2' */
    if (mainV03_56_B.Sum_mf > mainV03_56_P.Saturation2_UpperSat_d) {
      mainV03_56_B.Saturation2 = mainV03_56_P.Saturation2_UpperSat_d;
    } else if (mainV03_56_B.Sum_mf < mainV03_56_P.Saturation2_LowerSat_d) {
      mainV03_56_B.Saturation2 = mainV03_56_P.Saturation2_LowerSat_d;
    } else {
      mainV03_56_B.Saturation2 = mainV03_56_B.Sum_mf;
    }

    /* End of Saturate: '<S166>/Saturation2' */

    /* Sum: '<S166>/Sum3' */
    mainV03_56_B.Sum3_b = mainV03_56_B.Saturation2 - mainV03_56_B.ToBodyAxes7[1];

    /* Gain: '<S205>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_oc = mainV03_56_P.PIDController1_P_c *
      mainV03_56_B.Sum3_b;

    /* Integrator: '<S205>/Integrator' */
    mainV03_56_B.Integrator_lxr = mainV03_56_X.Integrator_CSTATE_m3;

    /* Gain: '<S205>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_m = mainV03_56_P.PIDController1_D_m *
      mainV03_56_B.Sum3_b;

    /* Integrator: '<S205>/Filter' */
    mainV03_56_B.Filter_gu = mainV03_56_X.Filter_CSTATE_kdq;

    /* Sum: '<S205>/SumD' */
    mainV03_56_B.SumD_hj = mainV03_56_B.DerivativeGain_m -
      mainV03_56_B.Filter_gu;

    /* Gain: '<S205>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_i = mainV03_56_P.PIDController1_N_c *
      mainV03_56_B.SumD_hj;

    /* Sum: '<S205>/Sum' */
    mainV03_56_B.Sum_lh = (mainV03_56_B.ProportionalGain_oc +
      mainV03_56_B.Integrator_lxr) + mainV03_56_B.FilterCoefficient_i;

    /* Saturate: '<S166>/Saturation1' */
    if (mainV03_56_B.Sum_lh > mainV03_56_P.Saturation1_UpperSat_o) {
      mainV03_56_B.Saturation1 = mainV03_56_P.Saturation1_UpperSat_o;
    } else if (mainV03_56_B.Sum_lh < mainV03_56_P.Saturation1_LowerSat_p) {
      mainV03_56_B.Saturation1 = mainV03_56_P.Saturation1_LowerSat_p;
    } else {
      mainV03_56_B.Saturation1 = mainV03_56_B.Sum_lh;
    }

    /* End of Saturate: '<S166>/Saturation1' */

    /* Sum: '<S166>/Sum6' */
    mainV03_56_B.Sum6_n = mainV03_56_B.Saturation1 -
      mainV03_56_B.Sensors.EulerMeas[0];

    /* Gain: '<S219>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_h = mainV03_56_P.ProportionalRoll_P_g *
      mainV03_56_B.Sum6_n;

    /* Integrator: '<S219>/Integrator' */
    mainV03_56_B.Integrator_p = mainV03_56_X.Integrator_CSTATE_gq;

    /* Gain: '<S219>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_k = mainV03_56_P.ProportionalRoll_D_e *
      mainV03_56_B.Sum6_n;

    /* Integrator: '<S219>/Filter' */
    mainV03_56_B.Filter_eg = mainV03_56_X.Filter_CSTATE_az;

    /* Sum: '<S219>/SumD' */
    mainV03_56_B.SumD_fj = mainV03_56_B.DerivativeGain_k -
      mainV03_56_B.Filter_eg;

    /* Gain: '<S219>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_mf = mainV03_56_P.ProportionalRoll_N_k *
      mainV03_56_B.SumD_fj;

    /* Sum: '<S219>/Sum' */
    mainV03_56_B.Sum_lq = (mainV03_56_B.ProportionalGain_h +
      mainV03_56_B.Integrator_p) + mainV03_56_B.FilterCoefficient_mf;

    /* Sum: '<S166>/Sum7' */
    mainV03_56_B.Sum7_b = mainV03_56_B.Sum_lq -
      mainV03_56_B.Sensors.OmegaMeas_body.pMeas;

    /* Gain: '<S215>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_le = mainV03_56_P.PIDRollRate_P_j *
      mainV03_56_B.Sum7_b;

    /* Integrator: '<S215>/Integrator' */
    mainV03_56_B.Integrator_gq = mainV03_56_X.Integrator_CSTATE_nb;

    /* Gain: '<S215>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_p = mainV03_56_P.PIDRollRate_D_l *
      mainV03_56_B.Sum7_b;

    /* Integrator: '<S215>/Filter' */
    mainV03_56_B.Filter_ia = mainV03_56_X.Filter_CSTATE_f;

    /* Sum: '<S215>/SumD' */
    mainV03_56_B.SumD_n = mainV03_56_B.DerivativeGain_p - mainV03_56_B.Filter_ia;

    /* Gain: '<S215>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_ps = mainV03_56_P.PIDRollRate_N_d *
      mainV03_56_B.SumD_n;

    /* Sum: '<S215>/Sum' */
    mainV03_56_B.Sum_je = (mainV03_56_B.ProportionalGain_le +
      mainV03_56_B.Integrator_gq) + mainV03_56_B.FilterCoefficient_ps;

    /* Gain: '<S225>/dRoll' */
    mainV03_56_B.dRoll = mainV03_56_P.dRoll_Gain_g * mainV03_56_B.Sum_je;

    /* Gain: '<S225>/Gain' */
    mainV03_56_B.Gain_p = mainV03_56_P.Gain_Gain_ci * mainV03_56_B.dRoll;

    /* Gain: '<S223>/dThrottle' */
    mainV03_56_B.dThrottle = mainV03_56_P.dThrottle_Gain_l *
      mainV03_56_B.ManualSwitch_o;

    /* Bias: '<S223>/Bias' */
    mainV03_56_B.Bias_g = mainV03_56_B.dThrottle + mainV03_56_P.Bias_Bias_o;

    /* ManualSwitch: '<S166>/Manual Switch6' incorporates:
     *  Constant: '<S166>/Constant4'
     */
    if (mainV03_56_P.ManualSwitch6_CurrentSetting_i == 1) {
      mainV03_56_B.ManualSwitch6 = mainV03_56_P.Constant4_Value_b;
    } else {
      mainV03_56_B.ManualSwitch6 = mainV03_56_B.ManualSwitch.yaw_cmd;
    }

    /* End of ManualSwitch: '<S166>/Manual Switch6' */

    /* Sum: '<S166>/Sum9' */
    mainV03_56_B.Sum9 = mainV03_56_B.ManualSwitch6 -
      mainV03_56_B.Sensors.EulerMeas[2];

    /* Gain: '<S220>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_ky = mainV03_56_P.ProportionalYaw1_P_p *
      mainV03_56_B.Sum9;

    /* Integrator: '<S220>/Integrator' */
    mainV03_56_B.Integrator_a = mainV03_56_X.Integrator_CSTATE_j;

    /* Gain: '<S220>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_b1 = mainV03_56_P.ProportionalYaw1_D_h *
      mainV03_56_B.Sum9;

    /* Integrator: '<S220>/Filter' */
    mainV03_56_B.Filter_er = mainV03_56_X.Filter_CSTATE_ok;

    /* Sum: '<S220>/SumD' */
    mainV03_56_B.SumD_ae = mainV03_56_B.DerivativeGain_b1 -
      mainV03_56_B.Filter_er;

    /* Gain: '<S220>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_ip = mainV03_56_P.ProportionalYaw1_N_j *
      mainV03_56_B.SumD_ae;

    /* Sum: '<S220>/Sum' */
    mainV03_56_B.Sum_asv = (mainV03_56_B.ProportionalGain_ky +
      mainV03_56_B.Integrator_a) + mainV03_56_B.FilterCoefficient_ip;

    /* Sum: '<S166>/Sum10' */
    mainV03_56_B.Sum10 = mainV03_56_B.Sum_asv -
      mainV03_56_B.Sensors.OmegaMeas_body.rMeas;

    /* Gain: '<S216>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_bf = mainV03_56_P.PIDYawRate1_P_g *
      mainV03_56_B.Sum10;

    /* Integrator: '<S216>/Integrator' */
    mainV03_56_B.Integrator_bh = mainV03_56_X.Integrator_CSTATE_h;

    /* Gain: '<S216>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_n5 = mainV03_56_P.PIDYawRate1_D_e *
      mainV03_56_B.Sum10;

    /* Integrator: '<S216>/Filter' */
    mainV03_56_B.Filter_nz = mainV03_56_X.Filter_CSTATE_or;

    /* Sum: '<S216>/SumD' */
    mainV03_56_B.SumD_bp = mainV03_56_B.DerivativeGain_n5 -
      mainV03_56_B.Filter_nz;

    /* Gain: '<S216>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_pt = mainV03_56_P.PIDYawRate1_N_h *
      mainV03_56_B.SumD_bp;

    /* Sum: '<S216>/Sum' */
    mainV03_56_B.Sum_fu = (mainV03_56_B.ProportionalGain_bf +
      mainV03_56_B.Integrator_bh) + mainV03_56_B.FilterCoefficient_pt;

    /* Gain: '<S226>/dYaw' */
    mainV03_56_B.dYaw = mainV03_56_P.dYaw_Gain_p * mainV03_56_B.Sum_fu;

    /* Saturate: '<S166>/Saturation7' */
    if (mainV03_56_B.Sum_oz > mainV03_56_P.Saturation7_UpperSat_d) {
      mainV03_56_B.Saturation7 = mainV03_56_P.Saturation7_UpperSat_d;
    } else if (mainV03_56_B.Sum_oz < mainV03_56_P.Saturation7_LowerSat_l) {
      mainV03_56_B.Saturation7 = mainV03_56_P.Saturation7_LowerSat_l;
    } else {
      mainV03_56_B.Saturation7 = mainV03_56_B.Sum_oz;
    }

    /* End of Saturate: '<S166>/Saturation7' */

    /* Sum: '<S166>/Sum15' */
    mainV03_56_B.Sum15 = mainV03_56_B.Saturation7 - mainV03_56_B.ToBodyAxes5[0];

    /* Gain: '<S204>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_lw = mainV03_56_P.PIDController_P_g *
      mainV03_56_B.Sum15;

    /* Integrator: '<S204>/Integrator' */
    mainV03_56_B.Integrator_av = mainV03_56_X.Integrator_CSTATE_b;

    /* Gain: '<S204>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_o = mainV03_56_P.PIDController_D_b *
      mainV03_56_B.Sum15;

    /* Integrator: '<S204>/Filter' */
    mainV03_56_B.Filter_n4 = mainV03_56_X.Filter_CSTATE_j;

    /* Sum: '<S204>/SumD' */
    mainV03_56_B.SumD_f5 = mainV03_56_B.DerivativeGain_o -
      mainV03_56_B.Filter_n4;

    /* Gain: '<S204>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_o5 = mainV03_56_P.PIDController_N_g *
      mainV03_56_B.SumD_f5;

    /* ManualSwitch: '<S166>/Manual Switch9' */
    if (mainV03_56_P.ManualSwitch9_CurrentSetting == 1) {
      /* ManualSwitch: '<S166>/Manual Switch11' incorporates:
       *  Constant: '<S166>/Constant1'
       */
      if (mainV03_56_P.ManualSwitch11_CurrentSetting == 1) {
        mainV03_56_B.ManualSwitch11 = mainV03_56_B.Sum19;
      } else {
        mainV03_56_B.ManualSwitch11 = mainV03_56_P.Constant1_Value_pn;
      }

      /* End of ManualSwitch: '<S166>/Manual Switch11' */
      mainV03_56_B.ManualSwitch9 = mainV03_56_B.ManualSwitch11;
    } else {
      /* Sum: '<S204>/Sum' */
      mainV03_56_B.Sum_iw = (mainV03_56_B.ProportionalGain_lw +
        mainV03_56_B.Integrator_av) + mainV03_56_B.FilterCoefficient_o5;

      /* Saturate: '<S166>/Saturation6' */
      if (mainV03_56_B.Sum_iw > mainV03_56_P.Saturation6_UpperSat_i) {
        mainV03_56_B.Saturation6 = mainV03_56_P.Saturation6_UpperSat_i;
      } else if (mainV03_56_B.Sum_iw < mainV03_56_P.Saturation6_LowerSat_f) {
        mainV03_56_B.Saturation6 = mainV03_56_P.Saturation6_LowerSat_f;
      } else {
        mainV03_56_B.Saturation6 = mainV03_56_B.Sum_iw;
      }

      /* End of Saturate: '<S166>/Saturation6' */
      mainV03_56_B.ManualSwitch9 = mainV03_56_B.Saturation6;
    }

    /* End of ManualSwitch: '<S166>/Manual Switch9' */

    /* Sum: '<S166>/Sum12' */
    mainV03_56_B.Sum12_m = mainV03_56_B.ManualSwitch9 -
      mainV03_56_B.Sensors.EulerMeas[1];

    /* Gain: '<S218>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_d = mainV03_56_P.ProportionalPitch2_P_i *
      mainV03_56_B.Sum12_m;

    /* Integrator: '<S218>/Integrator' */
    mainV03_56_B.Integrator_d2 = mainV03_56_X.Integrator_CSTATE_e5;

    /* Gain: '<S218>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_au = mainV03_56_P.ProportionalPitch2_D_f *
      mainV03_56_B.Sum12_m;

    /* Integrator: '<S218>/Filter' */
    mainV03_56_B.Filter_p = mainV03_56_X.Filter_CSTATE_mkr;

    /* Sum: '<S218>/SumD' */
    mainV03_56_B.SumD_e = mainV03_56_B.DerivativeGain_au - mainV03_56_B.Filter_p;

    /* Gain: '<S218>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_m4 = mainV03_56_P.ProportionalPitch2_N_b *
      mainV03_56_B.SumD_e;

    /* Sum: '<S218>/Sum' */
    mainV03_56_B.Sum_oy = (mainV03_56_B.ProportionalGain_d +
      mainV03_56_B.Integrator_d2) + mainV03_56_B.FilterCoefficient_m4;

    /* Sum: '<S166>/Sum13' */
    mainV03_56_B.PitchRateIn_g = mainV03_56_B.Sum_oy -
      mainV03_56_B.Sensors.OmegaMeas_body.qMeas;

    /* Gain: '<S214>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_at = mainV03_56_P.PIDPitchRate2_P_h *
      mainV03_56_B.PitchRateIn_g;

    /* Integrator: '<S214>/Integrator' */
    mainV03_56_B.Integrator_bf = mainV03_56_X.Integrator_CSTATE_ni;

    /* Gain: '<S214>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_hn = mainV03_56_P.PIDPitchRate2_D_d *
      mainV03_56_B.PitchRateIn_g;

    /* Integrator: '<S214>/Filter' */
    mainV03_56_B.Filter_o0m = mainV03_56_X.Filter_CSTATE_l;

    /* Sum: '<S214>/SumD' */
    mainV03_56_B.SumD_pe = mainV03_56_B.DerivativeGain_hn -
      mainV03_56_B.Filter_o0m;

    /* Gain: '<S214>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_dq = mainV03_56_P.PIDPitchRate2_N_j *
      mainV03_56_B.SumD_pe;

    /* Sum: '<S214>/Sum' */
    mainV03_56_B.Sum_mp = (mainV03_56_B.ProportionalGain_at +
      mainV03_56_B.Integrator_bf) + mainV03_56_B.FilterCoefficient_dq;

    /* Gain: '<S224>/dPitch' */
    mainV03_56_B.dPitch = mainV03_56_P.dPitch_Gain_g * mainV03_56_B.Sum_mp;

    /* Sum: '<S221>/Sum2' */
    mainV03_56_B.Throttle2noSaturation = ((mainV03_56_B.Gain_p +
      mainV03_56_B.Bias_g) + mainV03_56_B.dYaw) + mainV03_56_B.dPitch;

    /* Saturate: '<S221>/Saturation' */
    if (mainV03_56_B.Throttle2noSaturation > mainV03_56_P.Saturation_UpperSat_m)
    {
      mainV03_56_B.BusCreator1_g.Throttle2 = mainV03_56_P.Saturation_UpperSat_m;
    } else if (mainV03_56_B.Throttle2noSaturation <
               mainV03_56_P.Saturation_LowerSat_p) {
      mainV03_56_B.BusCreator1_g.Throttle2 = mainV03_56_P.Saturation_LowerSat_p;
    } else {
      mainV03_56_B.BusCreator1_g.Throttle2 = mainV03_56_B.Throttle2noSaturation;
    }

    /* End of Saturate: '<S221>/Saturation' */
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* Product: '<S224>/Divide' incorporates:
       *  Constant: '<S224>/Constant1'
       *  Constant: '<S224>/Constant2'
       *  Constant: '<S224>/motor2 position '
       *  Constant: '<S224>/motor2 position 1'
       *  Sum: '<S224>/Sum1'
       *  Sum: '<S224>/Sum2'
       */
      mainV03_56_B.Divide = (mainV03_56_P.motor2position_Value_o[0] -
        mainV03_56_P.Constant1_Value_bj[0]) /
        (mainV03_56_P.motor2position1_Value_n[0] -
         mainV03_56_P.Constant2_Value_c[0]);
    }

    /* Product: '<S224>/Product' */
    mainV03_56_B.Product_e4 = mainV03_56_B.dPitch * mainV03_56_B.Divide;

    /* Gain: '<S226>/Gain' */
    mainV03_56_B.Gain_g = mainV03_56_P.Gain_Gain_fs * mainV03_56_B.dYaw;

    /* Gain: '<S225>/Gain1' */
    mainV03_56_B.Gain1_n = mainV03_56_P.Gain1_Gain_o * mainV03_56_B.dRoll;

    /* Sum: '<S221>/Sum1' */
    mainV03_56_B.Throttle3noSaturation = ((mainV03_56_B.Product_e4 +
      mainV03_56_B.Bias_g) + mainV03_56_B.Gain_g) + mainV03_56_B.Gain1_n;

    /* Saturate: '<S221>/Saturation1' */
    if (mainV03_56_B.Throttle3noSaturation > mainV03_56_P.Saturation1_UpperSat_m)
    {
      mainV03_56_B.BusCreator1_g.Throttle3 = mainV03_56_P.Saturation1_UpperSat_m;
    } else if (mainV03_56_B.Throttle3noSaturation <
               mainV03_56_P.Saturation1_LowerSat_e) {
      mainV03_56_B.BusCreator1_g.Throttle3 = mainV03_56_P.Saturation1_LowerSat_e;
    } else {
      mainV03_56_B.BusCreator1_g.Throttle3 = mainV03_56_B.Throttle3noSaturation;
    }

    /* End of Saturate: '<S221>/Saturation1' */

    /* Product: '<S224>/Product1' */
    mainV03_56_B.Product1_b = mainV03_56_B.dPitch * mainV03_56_B.Divide;

    /* Sum: '<S221>/Sum3' */
    mainV03_56_B.Throttle4noSaturation = ((mainV03_56_B.dRoll +
      mainV03_56_B.Bias_g) + mainV03_56_B.dYaw) + mainV03_56_B.Product1_b;

    /* Saturate: '<S221>/Saturation2' */
    if (mainV03_56_B.Throttle4noSaturation > mainV03_56_P.Saturation2_UpperSat_i)
    {
      mainV03_56_B.BusCreator1_g.Throttle4 = mainV03_56_P.Saturation2_UpperSat_i;
    } else if (mainV03_56_B.Throttle4noSaturation <
               mainV03_56_P.Saturation2_LowerSat_p) {
      mainV03_56_B.BusCreator1_g.Throttle4 = mainV03_56_P.Saturation2_LowerSat_p;
    } else {
      mainV03_56_B.BusCreator1_g.Throttle4 = mainV03_56_B.Throttle4noSaturation;
    }

    /* End of Saturate: '<S221>/Saturation2' */

    /* Gain: '<S226>/Gain1' */
    mainV03_56_B.Gain1_b = mainV03_56_P.Gain1_Gain_f * mainV03_56_B.dYaw;

    /* Sum: '<S221>/Sum' */
    mainV03_56_B.Throttle5noSaturation = ((mainV03_56_B.dPitch +
      mainV03_56_B.Bias_g) + mainV03_56_B.Gain1_b) + mainV03_56_B.dRoll;

    /* Saturate: '<S221>/Saturation3' */
    if (mainV03_56_B.Throttle5noSaturation > mainV03_56_P.Saturation3_UpperSat_i)
    {
      mainV03_56_B.BusCreator1_g.Throttle5 = mainV03_56_P.Saturation3_UpperSat_i;
    } else if (mainV03_56_B.Throttle5noSaturation <
               mainV03_56_P.Saturation3_LowerSat_f) {
      mainV03_56_B.BusCreator1_g.Throttle5 = mainV03_56_P.Saturation3_LowerSat_f;
    } else {
      mainV03_56_B.BusCreator1_g.Throttle5 = mainV03_56_B.Throttle5noSaturation;
    }

    /* End of Saturate: '<S221>/Saturation3' */

    /* BusCreator: '<S221>/Bus Creator1' incorporates:
     *  Constant: '<S166>/Constant3'
     */
    mainV03_56_B.BusCreator1_g.Throttle1 = mainV03_56_P.Constant3_Value_n;
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* ManualSwitch: '<S166>/Manual Switch5' */
      if (mainV03_56_P.ManualSwitch5_CurrentSetting_i == 1) {
        mainV03_56_B.ManualSwitch5_l.Throttle1 = rtb_BusCreator1_Throttle1;
        mainV03_56_B.ManualSwitch5_l.Throttle2 = rtb_BusCreator1_Throttle2;
        mainV03_56_B.ManualSwitch5_l.Throttle3 = rtb_BusCreator1_Throttle3;
        mainV03_56_B.ManualSwitch5_l.Throttle4 = rtb_BusCreator1_Throttle4;
        mainV03_56_B.ManualSwitch5_l.Throttle5 = rtb_BusCreator1_Throttle5;
      } else {
        mainV03_56_B.ManualSwitch5_l.Throttle1 = rtb_BusCreator3_o_Throttle1;
        mainV03_56_B.ManualSwitch5_l.Throttle2 = rtb_BusCreator3_o_Throttle2;
        mainV03_56_B.ManualSwitch5_l.Throttle3 = rtb_BusCreator3_o_Throttle3;
        mainV03_56_B.ManualSwitch5_l.Throttle4 = rtb_BusCreator3_o_Throttle4;
        mainV03_56_B.ManualSwitch5_l.Throttle5 = rtb_BusCreator3_o_Throttle5;
      }

      /* End of ManualSwitch: '<S166>/Manual Switch5' */
    }

    /* ManualSwitch: '<S166>/Manual Switch3' */
    if (mainV03_56_P.ManualSwitch3_CurrentSetting_d == 1) {
      mainV03_56_B.ManualSwitch3 = mainV03_56_B.BusCreator1_g;
    } else {
      mainV03_56_B.ManualSwitch3 = mainV03_56_B.ManualSwitch5_l;
    }

    /* End of ManualSwitch: '<S166>/Manual Switch3' */

    /* Gain: '<S201>/Integral Gain' */
    mainV03_56_B.IntegralGain_m = mainV03_56_P.PHeight1_I *
      mainV03_56_B.Saturation11;

    /* Gain: '<S202>/Integral Gain' */
    mainV03_56_B.IntegralGain_e = mainV03_56_P.PIDAltitudeAcceleration_I_d *
      mainV03_56_B.Sum11;

    /* Sum: '<S166>/Sum' */
    mainV03_56_B.Sum_cs0 = mainV03_56_B.ManualSwitch.altitude_cmd -
      mainV03_56_B.Sensors.LLAMeas.AltitudeMeas_m;

    /* Gain: '<S217>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_au = mainV03_56_P.ProportionalAltitude_P_e *
      mainV03_56_B.Sum_cs0;

    /* Integrator: '<S217>/Integrator' */
    mainV03_56_B.Integrator_jk = mainV03_56_X.Integrator_CSTATE_c;

    /* Gain: '<S217>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_l = mainV03_56_P.ProportionalAltitude_D_h *
      mainV03_56_B.Sum_cs0;

    /* Integrator: '<S217>/Filter' */
    mainV03_56_B.Filter_m = mainV03_56_X.Filter_CSTATE_j2;

    /* Sum: '<S217>/SumD' */
    mainV03_56_B.SumD_a1 = mainV03_56_B.DerivativeGain_l - mainV03_56_B.Filter_m;

    /* Gain: '<S217>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_fx = mainV03_56_P.ProportionalAltitude_N_e *
      mainV03_56_B.SumD_a1;

    /* Sum: '<S217>/Sum' */
    mainV03_56_B.Sum_i = (mainV03_56_B.ProportionalGain_au +
                          mainV03_56_B.Integrator_jk) +
      mainV03_56_B.FilterCoefficient_fx;

    /* Saturate: '<S217>/Saturate' */
    if (mainV03_56_B.Sum_i >
        mainV03_56_P.ProportionalAltitude_UpperSaturationLimit_o) {
      mainV03_56_B.Saturate_g =
        mainV03_56_P.ProportionalAltitude_UpperSaturationLimit_o;
    } else if (mainV03_56_B.Sum_i <
               mainV03_56_P.ProportionalAltitude_LowerSaturationLimit_c) {
      mainV03_56_B.Saturate_g =
        mainV03_56_P.ProportionalAltitude_LowerSaturationLimit_c;
    } else {
      mainV03_56_B.Saturate_g = mainV03_56_B.Sum_i;
    }

    /* End of Saturate: '<S217>/Saturate' */

    /* Sum: '<S166>/Sum1' incorporates:
     *  UnaryMinus: '<S166>/Unary Minus'
     */
    mainV03_56_B.Sum1_hp = mainV03_56_B.Saturate_g -
      (-mainV03_56_B.Sensors.VearthMeas[2]);

    /* Gain: '<S203>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_bm = mainV03_56_P.PIDAltitudeRate_D_h *
      mainV03_56_B.Sum1_hp;

    /* Integrator: '<S203>/Filter' */
    mainV03_56_B.Filter_il = mainV03_56_X.Filter_CSTATE_br;

    /* Sum: '<S203>/SumD' */
    mainV03_56_B.SumD_ng = mainV03_56_B.DerivativeGain_bm -
      mainV03_56_B.Filter_il;

    /* Gain: '<S203>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_lp = mainV03_56_P.PIDAltitudeRate_N_k *
      mainV03_56_B.SumD_ng;

    /* Gain: '<S203>/Integral Gain' */
    mainV03_56_B.IntegralGain_fs = mainV03_56_P.PIDAltitudeRate_I_p *
      mainV03_56_B.Sum1_hp;

    /* Gain: '<S204>/Integral Gain' */
    mainV03_56_B.IntegralGain_o = mainV03_56_P.PIDController_I_b *
      mainV03_56_B.Sum15;

    /* Gain: '<S205>/Integral Gain' */
    mainV03_56_B.IntegralGain_g = mainV03_56_P.PIDController1_I_i *
      mainV03_56_B.Sum3_b;

    /* Gain: '<S206>/Integral Gain' */
    mainV03_56_B.IntegralGain_ay = mainV03_56_P.PIDController4_I_f *
      mainV03_56_B.PitchRateIn_j;

    /* Gain: '<S207>/Integral Gain' */
    mainV03_56_B.IntegralGain_p = mainV03_56_P.PIDController5_I_p *
      mainV03_56_B.Sum19;

    /* Gain: '<S208>/Integral Gain' */
    mainV03_56_B.IntegralGain_k = mainV03_56_P.PIDController6_I *
      mainV03_56_B.Sum21;

    /* Gain: '<S209>/Integral Gain' */
    mainV03_56_B.IntegralGain_ox = mainV03_56_P.PIDController7_I_a *
      mainV03_56_B.Sum22;

    /* Gain: '<S210>/Integral Gain' */
    mainV03_56_B.IntegralGain_gz = mainV03_56_P.PIDController8_I *
      mainV03_56_B.Sum23;

    /* Gain: '<S211>/Integral Gain' */
    mainV03_56_B.IntegralGain_jv = mainV03_56_P.PIDController9_I *
      mainV03_56_B.Sum24;

    /* Gain: '<S212>/Integral Gain' */
    mainV03_56_B.IntegralGain_bz = mainV03_56_P.PIDHeightController1_I_n *
      mainV03_56_B.Saturation9;

    /* Gain: '<S213>/Integral Gain' */
    mainV03_56_B.IntegralGain_n = mainV03_56_P.PIDLateralPosition_I_n *
      mainV03_56_B.Sum8;

    /* Gain: '<S214>/Integral Gain' */
    mainV03_56_B.IntegralGain_d1 = mainV03_56_P.PIDPitchRate2_I_l *
      mainV03_56_B.PitchRateIn_g;

    /* Gain: '<S215>/Integral Gain' */
    mainV03_56_B.IntegralGain_ix = mainV03_56_P.PIDRollRate_I_n *
      mainV03_56_B.Sum7_b;

    /* Gain: '<S216>/Integral Gain' */
    mainV03_56_B.IntegralGain_kz = mainV03_56_P.PIDYawRate1_I_j *
      mainV03_56_B.Sum10;

    /* Gain: '<S217>/Integral Gain' */
    mainV03_56_B.IntegralGain_b1 = mainV03_56_P.ProportionalAltitude_I_h *
      mainV03_56_B.Sum_cs0;

    /* Gain: '<S218>/Integral Gain' */
    mainV03_56_B.IntegralGain_da = mainV03_56_P.ProportionalPitch2_I_d *
      mainV03_56_B.Sum12_m;

    /* Gain: '<S219>/Integral Gain' */
    mainV03_56_B.IntegralGain_eu = mainV03_56_P.ProportionalRoll_I_b *
      mainV03_56_B.Sum6_n;

    /* Gain: '<S220>/Integral Gain' */
    mainV03_56_B.IntegralGain_ao = mainV03_56_P.ProportionalYaw1_I_a *
      mainV03_56_B.Sum9;
  }

  /* End of Outputs for SubSystem: '<S3>/Quadcopter --> Fixed-Wing' */

  /* Product: '<S3>/Product5' */
  mainV03_56_B.Product5[0] = (real_T)mainV03_56_B.Compare_j *
    mainV03_56_B.BusCreator2.deltae;
  mainV03_56_B.Product5[1] = (real_T)mainV03_56_B.Compare_j *
    mainV03_56_B.BusCreator2.deltar;
  mainV03_56_B.Product5[2] = (real_T)mainV03_56_B.Compare_j *
    mainV03_56_B.BusCreator2.deltafr;
  mainV03_56_B.Product5[3] = (real_T)mainV03_56_B.Compare_j *
    mainV03_56_B.BusCreator2.deltafl;

  /* RelationalOperator: '<S230>/Compare' incorporates:
   *  Constant: '<S230>/Constant'
   */
  mainV03_56_B.Compare_g = (mainV03_56_B.DataTypeConversion ==
    mainV03_56_P.CompareToConstant2_const);
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* SignalConversion: '<S3>/HiddenBuf_InsertedFor_Fixed-Wing Climb_at_inport_2' */
    mainV03_56_B.HiddenBuf_InsertedFor_FixedWingClimb_at_inport_2 =
      mainV03_56_B.Compare_g;

    /* Outputs for Enabled SubSystem: '<S3>/Fixed-Wing Climb' incorporates:
     *  EnablePort: '<S164>/Enable'
     */
    if (rtmIsMajorTimeStep(mainV03_56_M)) {
      if (mainV03_56_B.HiddenBuf_InsertedFor_FixedWingClimb_at_inport_2) {
        if (!mainV03_56_DW.FixedWingClimb_MODE) {
          mainV03_56_DW.FixedWingClimb_MODE = true;
        }
      } else {
        if (mainV03_56_DW.FixedWingClimb_MODE) {
          mainV03_56_DW.FixedWingClimb_MODE = false;
        }
      }
    }

    /* End of Outputs for SubSystem: '<S3>/Fixed-Wing Climb' */
  }

  /* Outputs for Enabled SubSystem: '<S3>/Fixed-Wing Climb' incorporates:
   *  EnablePort: '<S164>/Enable'
   */
  if (mainV03_56_DW.FixedWingClimb_MODE) {
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* Gain: '<S164>/Gain' incorporates:
       *  Constant: '<S164>/Constant9'
       */
      mainV03_56_B.actuators_f.deltae = mainV03_56_P.Gain_Gain_b *
        mainV03_56_P.deltae_degrees;

      /* Gain: '<S164>/Gain1' incorporates:
       *  Constant: '<S164>/Constant10'
       */
      mainV03_56_B.actuators_f.deltar = mainV03_56_P.Gain1_Gain_p *
        mainV03_56_P.deltar_degrees;

      /* Gain: '<S164>/Gain2' incorporates:
       *  Constant: '<S164>/Constant1'
       */
      mainV03_56_B.actuators_f.deltafr = mainV03_56_P.Gain2_Gain_m *
        mainV03_56_P.deltafr_degrees;

      /* Gain: '<S164>/Gain3' incorporates:
       *  Constant: '<S164>/Constant2'
       */
      mainV03_56_B.actuators_f.deltafl = mainV03_56_P.Gain3_Gain_d *
        mainV03_56_P.deltafl_degrees;
    }

    /* DotProduct: '<S164>/Dot Product' */
    riseValLimit = 0.0;
    for (rtb_Sum1_eh = 0; rtb_Sum1_eh < 3; rtb_Sum1_eh++) {
      /* Product: '<S164>/To Body Axes1' */
      mainV03_56_B.ToBodyAxes1[rtb_Sum1_eh] = 0.0;
      mainV03_56_B.ToBodyAxes1[rtb_Sum1_eh] +=
        mainV03_56_B.Sensors.DCMMeas_body_earth[rtb_Sum1_eh] *
        mainV03_56_B.Sensors.VearthMeas[0];
      mainV03_56_B.ToBodyAxes1[rtb_Sum1_eh] +=
        mainV03_56_B.Sensors.DCMMeas_body_earth[rtb_Sum1_eh + 3] *
        mainV03_56_B.Sensors.VearthMeas[1];
      mainV03_56_B.ToBodyAxes1[rtb_Sum1_eh] +=
        mainV03_56_B.Sensors.DCMMeas_body_earth[rtb_Sum1_eh + 6] *
        mainV03_56_B.Sensors.VearthMeas[2];

      /* DotProduct: '<S164>/Dot Product' */
      riseValLimit += mainV03_56_B.ToBodyAxes1[rtb_Sum1_eh] *
        mainV03_56_B.ToBodyAxes1[rtb_Sum1_eh];
    }

    /* Sum: '<S164>/Sum7' incorporates:
     *  Constant: '<S164>/Constant18'
     *  DotProduct: '<S164>/Dot Product'
     *  Sqrt: '<S164>/Sqrt'
     */
    mainV03_56_B.VbodyIn = mainV03_56_P.Constant18_Value - sqrt(riseValLimit);

    /* Gain: '<S176>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_a3 = mainV03_56_P.PIDController7_P *
      mainV03_56_B.VbodyIn;

    /* Integrator: '<S176>/Integrator' */
    mainV03_56_B.Integrator_ms = mainV03_56_X.Integrator_CSTATE_dg;

    /* Gain: '<S176>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_nt = mainV03_56_P.PIDController7_D *
      mainV03_56_B.VbodyIn;

    /* Integrator: '<S176>/Filter' */
    mainV03_56_B.Filter_mm = mainV03_56_X.Filter_CSTATE_k1;

    /* Sum: '<S176>/SumD' */
    mainV03_56_B.SumD_gx = mainV03_56_B.DerivativeGain_nt -
      mainV03_56_B.Filter_mm;

    /* Gain: '<S176>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_ge = mainV03_56_P.PIDController7_N *
      mainV03_56_B.SumD_gx;

    /* Sum: '<S176>/Sum' */
    mainV03_56_B.Sum_bi = (mainV03_56_B.ProportionalGain_a3 +
      mainV03_56_B.Integrator_ms) + mainV03_56_B.FilterCoefficient_ge;

    /* Saturate: '<S176>/Saturate' */
    if (mainV03_56_B.Sum_bi > mainV03_56_P.PIDController7_UpperSaturationLimit)
    {
      mainV03_56_B.Saturate_j = mainV03_56_P.PIDController7_UpperSaturationLimit;
    } else if (mainV03_56_B.Sum_bi <
               mainV03_56_P.PIDController7_LowerSaturationLimit) {
      mainV03_56_B.Saturate_j = mainV03_56_P.PIDController7_LowerSaturationLimit;
    } else {
      mainV03_56_B.Saturate_j = mainV03_56_B.Sum_bi;
    }

    /* End of Saturate: '<S176>/Saturate' */

    /* Sum: '<S164>/Sum1' incorporates:
     *  Constant: '<S164>/Constant15'
     */
    mainV03_56_B.Sum1_k = mainV03_56_P.Constant15_Value -
      mainV03_56_B.Sensors.LLAMeas.AltitudeMeas_m;

    /* Saturate: '<S164>/Saturation' */
    if (mainV03_56_B.Sum1_k > mainV03_56_P.Saturation_UpperSat) {
      mainV03_56_B.Saturation_m = mainV03_56_P.Saturation_UpperSat;
    } else if (mainV03_56_B.Sum1_k < mainV03_56_P.Saturation_LowerSat) {
      mainV03_56_B.Saturation_m = mainV03_56_P.Saturation_LowerSat;
    } else {
      mainV03_56_B.Saturation_m = mainV03_56_B.Sum1_k;
    }

    /* End of Saturate: '<S164>/Saturation' */

    /* Gain: '<S169>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_c1 = mainV03_56_P.PHeight_P *
      mainV03_56_B.Saturation_m;

    /* Integrator: '<S169>/Integrator' */
    mainV03_56_B.Integrator_la = mainV03_56_X.Integrator_CSTATE_ke;

    /* Gain: '<S169>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_kp = mainV03_56_P.PHeight_D *
      mainV03_56_B.Saturation_m;

    /* Integrator: '<S169>/Filter' */
    mainV03_56_B.Filter_m2 = mainV03_56_X.Filter_CSTATE_ba;

    /* Sum: '<S169>/SumD' */
    mainV03_56_B.SumD_me = mainV03_56_B.DerivativeGain_kp -
      mainV03_56_B.Filter_m2;

    /* Gain: '<S169>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_f0 = mainV03_56_P.PHeight_N *
      mainV03_56_B.SumD_me;

    /* Sum: '<S169>/Sum' */
    mainV03_56_B.Sum_gg = (mainV03_56_B.ProportionalGain_c1 +
      mainV03_56_B.Integrator_la) + mainV03_56_B.FilterCoefficient_f0;

    /* Saturate: '<S169>/Saturate' */
    if (mainV03_56_B.Sum_gg > mainV03_56_P.PHeight_UpperSaturationLimit) {
      mainV03_56_B.Saturate_o = mainV03_56_P.PHeight_UpperSaturationLimit;
    } else if (mainV03_56_B.Sum_gg < mainV03_56_P.PHeight_LowerSaturationLimit)
    {
      mainV03_56_B.Saturate_o = mainV03_56_P.PHeight_LowerSaturationLimit;
    } else {
      mainV03_56_B.Saturate_o = mainV03_56_B.Sum_gg;
    }

    /* End of Saturate: '<S169>/Saturate' */

    /* ManualSwitch: '<S164>/Manual Switch' */
    if (mainV03_56_P.ManualSwitch_CurrentSetting == 1) {
      mainV03_56_B.ManualSwitch_d = 0.0;
    } else {
      mainV03_56_B.ManualSwitch_d = mainV03_56_B.Saturate_o;
    }

    /* End of ManualSwitch: '<S164>/Manual Switch' */

    /* Sum: '<S164>/Sum12' */
    mainV03_56_B.Sum12_o = mainV03_56_B.ManualSwitch_d -
      mainV03_56_B.Sensors.alphaMeas;

    /* Gain: '<S171>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_cd = mainV03_56_P.PIDController1_P *
      mainV03_56_B.Sum12_o;

    /* Integrator: '<S171>/Integrator' */
    mainV03_56_B.Integrator_ew = mainV03_56_X.Integrator_CSTATE_ih;

    /* Gain: '<S171>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_fn = mainV03_56_P.PIDController1_D *
      mainV03_56_B.Sum12_o;

    /* Integrator: '<S171>/Filter' */
    mainV03_56_B.Filter_b4 = mainV03_56_X.Filter_CSTATE_lp;

    /* Sum: '<S171>/SumD' */
    mainV03_56_B.SumD_fjm = mainV03_56_B.DerivativeGain_fn -
      mainV03_56_B.Filter_b4;

    /* Gain: '<S171>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_n0 = mainV03_56_P.PIDController1_N *
      mainV03_56_B.SumD_fjm;

    /* Sum: '<S171>/Sum' */
    mainV03_56_B.Sum_j5 = (mainV03_56_B.ProportionalGain_cd +
      mainV03_56_B.Integrator_ew) + mainV03_56_B.FilterCoefficient_n0;

    /* Sum: '<S164>/Sum13' */
    mainV03_56_B.PitchRateIn_h = mainV03_56_B.Sum_j5 -
      mainV03_56_B.Sensors.OmegaMeas_body.qMeas;

    /* Gain: '<S170>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_j = mainV03_56_P.PIDController_P *
      mainV03_56_B.PitchRateIn_h;

    /* Integrator: '<S170>/Integrator' */
    mainV03_56_B.Integrator_lp4 = mainV03_56_X.Integrator_CSTATE_nm;

    /* Gain: '<S170>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_m4 = mainV03_56_P.PIDController_D *
      mainV03_56_B.PitchRateIn_h;

    /* Integrator: '<S170>/Filter' */
    mainV03_56_B.Filter_fn = mainV03_56_X.Filter_CSTATE_dt;

    /* Sum: '<S170>/SumD' */
    mainV03_56_B.SumD_na = mainV03_56_B.DerivativeGain_m4 -
      mainV03_56_B.Filter_fn;

    /* Gain: '<S170>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_jh = mainV03_56_P.PIDController_N *
      mainV03_56_B.SumD_na;

    /* Sum: '<S170>/Sum' */
    mainV03_56_B.Sum_it = (mainV03_56_B.ProportionalGain_j +
      mainV03_56_B.Integrator_lp4) + mainV03_56_B.FilterCoefficient_jh;

    /* Saturate: '<S170>/Saturate' */
    if (mainV03_56_B.Sum_it > mainV03_56_P.PIDController_UpperSaturationLimit) {
      mainV03_56_B.Saturate_ia = mainV03_56_P.PIDController_UpperSaturationLimit;
    } else if (mainV03_56_B.Sum_it <
               mainV03_56_P.PIDController_LowerSaturationLimit) {
      mainV03_56_B.Saturate_ia = mainV03_56_P.PIDController_LowerSaturationLimit;
    } else {
      mainV03_56_B.Saturate_ia = mainV03_56_B.Sum_it;
    }

    /* End of Saturate: '<S170>/Saturate' */
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* Gain: '<S164>/Gain4' incorporates:
       *  Constant: '<S164>/Constant'
       */
      mainV03_56_B.deltae_f = mainV03_56_P.Gain4_Gain * mainV03_56_P.delta_e;
    }

    /* ManualSwitch: '<S164>/Manual Switch2' */
    if (mainV03_56_P.ManualSwitch2_CurrentSetting == 1) {
      /* ManualSwitch: '<S164>/Manual Switch7' */
      if (mainV03_56_P.ManualSwitch7_CurrentSetting == 1) {
        mainV03_56_B.ManualSwitch7 = mainV03_56_B.Saturate_j;
      } else {
        mainV03_56_B.ManualSwitch7 = mainV03_56_B.Saturate_ia;
      }

      /* End of ManualSwitch: '<S164>/Manual Switch7' */
      mainV03_56_B.BusCreator1_c.deltae = mainV03_56_B.ManualSwitch7;
    } else {
      mainV03_56_B.BusCreator1_c.deltae = mainV03_56_B.deltae_f;
    }

    /* End of ManualSwitch: '<S164>/Manual Switch2' */

    /* Sum: '<S164>/Sum3' incorporates:
     *  Constant: '<S164>/Constant14'
     */
    mainV03_56_B.Sum3_i = mainV03_56_P.Constant14_Value -
      mainV03_56_B.Sensors.EulerMeas[2];

    /* Gain: '<S173>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_i2 = mainV03_56_P.PIDController3_P *
      mainV03_56_B.Sum3_i;

    /* Integrator: '<S173>/Integrator' */
    mainV03_56_B.Integrator_et = mainV03_56_X.Integrator_CSTATE_og;

    /* Gain: '<S173>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_me = mainV03_56_P.PIDController3_D *
      mainV03_56_B.Sum3_i;

    /* Integrator: '<S173>/Filter' */
    mainV03_56_B.Filter_p4 = mainV03_56_X.Filter_CSTATE_iv;

    /* Sum: '<S173>/SumD' */
    mainV03_56_B.SumD_hb = mainV03_56_B.DerivativeGain_me -
      mainV03_56_B.Filter_p4;

    /* Gain: '<S173>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_d3 = mainV03_56_P.PIDController3_N *
      mainV03_56_B.SumD_hb;

    /* Sum: '<S173>/Sum' */
    mainV03_56_B.Sum_i1 = (mainV03_56_B.ProportionalGain_i2 +
      mainV03_56_B.Integrator_et) + mainV03_56_B.FilterCoefficient_d3;

    /* Sum: '<S164>/Sum2' */
    mainV03_56_B.Sum2_aq = mainV03_56_B.Sum_i1 -
      mainV03_56_B.Sensors.OmegaMeas_body.rMeas;

    /* Gain: '<S172>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_ce = mainV03_56_P.PIDController2_P *
      mainV03_56_B.Sum2_aq;

    /* Integrator: '<S172>/Integrator' */
    mainV03_56_B.Integrator_o = mainV03_56_X.Integrator_CSTATE_ht;

    /* Gain: '<S172>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_jw = mainV03_56_P.PIDController2_D *
      mainV03_56_B.Sum2_aq;

    /* Integrator: '<S172>/Filter' */
    mainV03_56_B.Filter_c = mainV03_56_X.Filter_CSTATE_k5;

    /* Sum: '<S172>/SumD' */
    mainV03_56_B.SumD_jm = mainV03_56_B.DerivativeGain_jw -
      mainV03_56_B.Filter_c;

    /* Gain: '<S172>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_hr = mainV03_56_P.PIDController2_N *
      mainV03_56_B.SumD_jm;
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* Gain: '<S164>/Gain5' incorporates:
       *  Constant: '<S164>/Constant11'
       */
      mainV03_56_B.deltar_k = mainV03_56_P.Gain5_Gain *
        mainV03_56_P.Constant11_Value;
    }

    /* ManualSwitch: '<S164>/Manual Switch3' */
    if (mainV03_56_P.ManualSwitch3_CurrentSetting == 1) {
      /* Sum: '<S172>/Sum' */
      mainV03_56_B.Sum_ggh = (mainV03_56_B.ProportionalGain_ce +
        mainV03_56_B.Integrator_o) + mainV03_56_B.FilterCoefficient_hr;
      mainV03_56_B.BusCreator1_c.deltar = mainV03_56_B.Sum_ggh;
    } else {
      mainV03_56_B.BusCreator1_c.deltar = mainV03_56_B.deltar_k;
    }

    /* End of ManualSwitch: '<S164>/Manual Switch3' */

    /* Sum: '<S164>/Sum5' incorporates:
     *  Constant: '<S164>/Constant16'
     */
    mainV03_56_B.Sum5_c = mainV03_56_P.Constant16_Value -
      mainV03_56_B.Sensors.EulerMeas[0];

    /* Gain: '<S175>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_bl = mainV03_56_P.PIDController5_P *
      mainV03_56_B.Sum5_c;

    /* Integrator: '<S175>/Integrator' */
    mainV03_56_B.Integrator_lr = mainV03_56_X.Integrator_CSTATE_cw;

    /* Gain: '<S175>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_hj = mainV03_56_P.PIDController5_D *
      mainV03_56_B.Sum5_c;

    /* Integrator: '<S175>/Filter' */
    mainV03_56_B.Filter_dp = mainV03_56_X.Filter_CSTATE_lo;

    /* Sum: '<S175>/SumD' */
    mainV03_56_B.SumD_kp = mainV03_56_B.DerivativeGain_hj -
      mainV03_56_B.Filter_dp;

    /* Gain: '<S175>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_ok = mainV03_56_P.PIDController5_N *
      mainV03_56_B.SumD_kp;

    /* Sum: '<S175>/Sum' */
    mainV03_56_B.Sum_mm = (mainV03_56_B.ProportionalGain_bl +
      mainV03_56_B.Integrator_lr) + mainV03_56_B.FilterCoefficient_ok;

    /* Sum: '<S164>/Sum4' */
    mainV03_56_B.Sum4_k = mainV03_56_B.Sum_mm -
      mainV03_56_B.Sensors.OmegaMeas_body.pMeas;

    /* Gain: '<S174>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_j2 = mainV03_56_P.PIDController4_P *
      mainV03_56_B.Sum4_k;

    /* Integrator: '<S174>/Integrator' */
    mainV03_56_B.Integrator_ey = mainV03_56_X.Integrator_CSTATE_ki;

    /* Gain: '<S174>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_ic = mainV03_56_P.PIDController4_D *
      mainV03_56_B.Sum4_k;

    /* Integrator: '<S174>/Filter' */
    mainV03_56_B.Filter_eo = mainV03_56_X.Filter_CSTATE_kc;

    /* Sum: '<S174>/SumD' */
    mainV03_56_B.SumD_pf = mainV03_56_B.DerivativeGain_ic -
      mainV03_56_B.Filter_eo;

    /* Gain: '<S174>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_cv = mainV03_56_P.PIDController4_N *
      mainV03_56_B.SumD_pf;

    /* Sum: '<S174>/Sum' */
    mainV03_56_B.Sum_gf = (mainV03_56_B.ProportionalGain_j2 +
      mainV03_56_B.Integrator_ey) + mainV03_56_B.FilterCoefficient_cv;
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* Gain: '<S164>/Gain6' incorporates:
       *  Constant: '<S164>/Constant12'
       */
      mainV03_56_B.deltafr_c = mainV03_56_P.Gain6_Gain *
        mainV03_56_P.Constant12_Value;

      /* Gain: '<S164>/Gain7' incorporates:
       *  Constant: '<S164>/Constant13'
       */
      mainV03_56_B.deltafl_p = mainV03_56_P.Gain7_Gain *
        mainV03_56_P.Constant13_Value;
    }

    /* ManualSwitch: '<S164>/Manual Switch4' */
    if (mainV03_56_P.ManualSwitch4_CurrentSetting == 1) {
      /* Sum: '<S164>/Sum6' incorporates:
       *  Constant: '<S164>/Constant17'
       */
      mainV03_56_B.Sum6_g = mainV03_56_B.Sum_gf + mainV03_56_P.Constant17_Value;
      mainV03_56_B.BusCreator1_c.deltafr = mainV03_56_B.Sum6_g;
    } else {
      mainV03_56_B.BusCreator1_c.deltafr = mainV03_56_B.deltafr_c;
    }

    /* End of ManualSwitch: '<S164>/Manual Switch4' */

    /* ManualSwitch: '<S164>/Manual Switch6' */
    if (mainV03_56_P.ManualSwitch6_CurrentSetting == 1) {
      /* Gain: '<S164>/Gain8' */
      mainV03_56_B.Gain8_i = mainV03_56_P.Gain8_Gain * mainV03_56_B.Sum_gf;

      /* Sum: '<S164>/Sum' incorporates:
       *  Constant: '<S164>/Constant17'
       */
      mainV03_56_B.Sum_ast = mainV03_56_B.Gain8_i +
        mainV03_56_P.Constant17_Value;
      mainV03_56_B.BusCreator1_c.deltafl = mainV03_56_B.Sum_ast;
    } else {
      mainV03_56_B.BusCreator1_c.deltafl = mainV03_56_B.deltafl_p;
    }

    /* End of ManualSwitch: '<S164>/Manual Switch6' */

    /* ManualSwitch: '<S164>/Manual Switch1' */
    if (mainV03_56_P.ManualSwitch1_CurrentSetting == 1) {
      mainV03_56_B.ManualSwitch1_n = mainV03_56_B.BusCreator1_c;
    } else {
      mainV03_56_B.ManualSwitch1_n = mainV03_56_B.actuators_f;
    }

    /* End of ManualSwitch: '<S164>/Manual Switch1' */

    /* Gain: '<S169>/Integral Gain' */
    mainV03_56_B.IntegralGain_m3 = mainV03_56_P.PHeight_I *
      mainV03_56_B.Saturation_m;

    /* Gain: '<S170>/Integral Gain' */
    mainV03_56_B.IntegralGain_jo = mainV03_56_P.PIDController_I *
      mainV03_56_B.PitchRateIn_h;

    /* Gain: '<S171>/Integral Gain' */
    mainV03_56_B.IntegralGain_ah = mainV03_56_P.PIDController1_I *
      mainV03_56_B.Sum12_o;

    /* Gain: '<S172>/Integral Gain' */
    mainV03_56_B.IntegralGain_fk = mainV03_56_P.PIDController2_I *
      mainV03_56_B.Sum2_aq;

    /* Gain: '<S173>/Integral Gain' */
    mainV03_56_B.IntegralGain_pd = mainV03_56_P.PIDController3_I *
      mainV03_56_B.Sum3_i;

    /* Gain: '<S174>/Integral Gain' */
    mainV03_56_B.IntegralGain_hb = mainV03_56_P.PIDController4_I *
      mainV03_56_B.Sum4_k;

    /* Gain: '<S175>/Integral Gain' */
    mainV03_56_B.IntegralGain_a3 = mainV03_56_P.PIDController5_I *
      mainV03_56_B.Sum5_c;

    /* Gain: '<S176>/Integral Gain' */
    mainV03_56_B.IntegralGain_cv = mainV03_56_P.PIDController7_I *
      mainV03_56_B.VbodyIn;
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* ManualSwitch: '<S164>/Manual Switch5' incorporates:
       *  BusCreator: '<S164>/Bus Creator3'
       *  Constant: '<S164>/Constant3'
       *  Constant: '<S164>/Constant4'
       *  Constant: '<S164>/Constant5'
       *  Constant: '<S164>/Constant7'
       *  Constant: '<S164>/Constant8'
       *  Constant: '<S164>/deltat_10'
       *  Constant: '<S164>/deltat_6'
       *  Constant: '<S164>/deltat_7'
       *  Constant: '<S164>/deltat_8'
       *  Constant: '<S164>/deltat_9'
       */
      if (mainV03_56_P.ManualSwitch5_CurrentSetting == 1) {
        mainV03_56_B.ManualSwitch5_j.Throttle1 = mainV03_56_P.Constant4_Value;
        mainV03_56_B.ManualSwitch5_j.Throttle2 = mainV03_56_P.Constant5_Value;
        mainV03_56_B.ManualSwitch5_j.Throttle3 = mainV03_56_P.Constant7_Value;
        mainV03_56_B.ManualSwitch5_j.Throttle4 = mainV03_56_P.Constant8_Value;
        mainV03_56_B.ManualSwitch5_j.Throttle5 = mainV03_56_P.Constant3_Value_m;
      } else {
        mainV03_56_B.ManualSwitch5_j.Throttle1 = mainV03_56_P.deltat_1;
        mainV03_56_B.ManualSwitch5_j.Throttle2 = mainV03_56_P.deltat_2;
        mainV03_56_B.ManualSwitch5_j.Throttle3 = mainV03_56_P.deltat_3;
        mainV03_56_B.ManualSwitch5_j.Throttle4 = mainV03_56_P.deltat_4;
        mainV03_56_B.ManualSwitch5_j.Throttle5 = mainV03_56_P.deltat_5;
      }

      /* End of ManualSwitch: '<S164>/Manual Switch5' */
    }
  }

  /* End of Outputs for SubSystem: '<S3>/Fixed-Wing Climb' */

  /* Product: '<S3>/Product3' */
  mainV03_56_B.Product3[0] = (real_T)mainV03_56_B.Compare_g *
    mainV03_56_B.ManualSwitch1_n.deltae;
  mainV03_56_B.Product3[1] = (real_T)mainV03_56_B.Compare_g *
    mainV03_56_B.ManualSwitch1_n.deltar;
  mainV03_56_B.Product3[2] = (real_T)mainV03_56_B.Compare_g *
    mainV03_56_B.ManualSwitch1_n.deltafr;
  mainV03_56_B.Product3[3] = (real_T)mainV03_56_B.Compare_g *
    mainV03_56_B.ManualSwitch1_n.deltafl;

  /* RelationalOperator: '<S231>/Compare' incorporates:
   *  Constant: '<S231>/Constant'
   */
  rtb_Compare_a = (mainV03_56_B.DataTypeConversion ==
                   mainV03_56_P.CompareToConstant3_const);

  /* Outputs for Enabled SubSystem: '<S3>/Fixed-Wing - Cruise' incorporates:
   *  EnablePort: '<S163>/Enable'
   */
  if (rtmIsMajorTimeStep(mainV03_56_M)) {
    if (rtb_Compare_a) {
      if (!mainV03_56_DW.FixedWingCruise_MODE) {
        mainV03_56_DW.FixedWingCruise_MODE = true;
      }
    } else {
      if (mainV03_56_DW.FixedWingCruise_MODE) {
        mainV03_56_DW.FixedWingCruise_MODE = false;
      }
    }
  }

  if (mainV03_56_DW.FixedWingCruise_MODE && (rtmIsMajorTimeStep(mainV03_56_M) &&
       mainV03_56_M->Timing.TaskCounters.TID[1] == 0)) {
    /* BusCreator: '<S163>/Bus Creator' incorporates:
     *  Constant: '<S163>/Constant1'
     *  Constant: '<S163>/Constant10'
     *  Constant: '<S163>/Constant2'
     *  Constant: '<S163>/Constant9'
     *  Gain: '<S163>/Gain'
     *  Gain: '<S163>/Gain1'
     *  Gain: '<S163>/Gain2'
     *  Gain: '<S163>/Gain3'
     */
    mainV03_56_B.Actuators.deltae = mainV03_56_P.Gain_Gain_i *
      mainV03_56_P.Constant9_Value;
    mainV03_56_B.Actuators.deltar = mainV03_56_P.Gain1_Gain_i *
      mainV03_56_P.Constant10_Value;
    mainV03_56_B.Actuators.deltafr = mainV03_56_P.Gain2_Gain_g *
      mainV03_56_P.Constant1_Value_b;
    mainV03_56_B.Actuators.deltafl = mainV03_56_P.Gain3_Gain *
      mainV03_56_P.Constant2_Value_h;

    /* InitialCondition: '<S163>/IC' incorporates:
     *  Constant: '<S163>/Constant21'
     */
    if (mainV03_56_DW.IC_FirstOutputTime_m) {
      mainV03_56_DW.IC_FirstOutputTime_m = false;
      mainV03_56_B.Throttle1 = mainV03_56_P.IC_Value;
    } else {
      mainV03_56_B.Throttle1 = mainV03_56_P.Constant21_Value;
    }

    /* End of InitialCondition: '<S163>/IC' */

    /* InitialCondition: '<S163>/IC1' incorporates:
     *  Constant: '<S163>/Constant22'
     */
    if (mainV03_56_DW.IC1_FirstOutputTime_c) {
      mainV03_56_DW.IC1_FirstOutputTime_c = false;
      mainV03_56_B.Throttle2_g = mainV03_56_P.IC1_Value;
    } else {
      mainV03_56_B.Throttle2_g = mainV03_56_P.Constant22_Value;
    }

    /* End of InitialCondition: '<S163>/IC1' */

    /* InitialCondition: '<S163>/IC2' incorporates:
     *  Constant: '<S163>/Constant23'
     */
    if (mainV03_56_DW.IC2_FirstOutputTime_do) {
      mainV03_56_DW.IC2_FirstOutputTime_do = false;
      mainV03_56_B.Throttle3_jz = mainV03_56_P.IC2_Value;
    } else {
      mainV03_56_B.Throttle3_jz = mainV03_56_P.Constant23_Value;
    }

    /* End of InitialCondition: '<S163>/IC2' */

    /* InitialCondition: '<S163>/IC3' incorporates:
     *  Constant: '<S163>/Constant24'
     */
    if (mainV03_56_DW.IC3_FirstOutputTime_m) {
      mainV03_56_DW.IC3_FirstOutputTime_m = false;
      mainV03_56_B.Throttle4_j = mainV03_56_P.IC3_Value;
    } else {
      mainV03_56_B.Throttle4_j = mainV03_56_P.Constant24_Value;
    }

    /* End of InitialCondition: '<S163>/IC3' */

    /* InitialCondition: '<S163>/IC4' incorporates:
     *  Constant: '<S163>/Constant20'
     */
    if (mainV03_56_DW.IC4_FirstOutputTime_c) {
      mainV03_56_DW.IC4_FirstOutputTime_c = false;
      mainV03_56_B.Throttle5_a = mainV03_56_P.IC4_Value;
    } else {
      mainV03_56_B.Throttle5_a = mainV03_56_P.Constant20_Value;
    }

    /* End of InitialCondition: '<S163>/IC4' */

    /* BusCreator: '<S163>/Bus Creator1' */
    mainV03_56_B.Throttle_m.Throttle1 = mainV03_56_B.Throttle1;
    mainV03_56_B.Throttle_m.Throttle2 = mainV03_56_B.Throttle2_g;
    mainV03_56_B.Throttle_m.Throttle3 = mainV03_56_B.Throttle3_jz;
    mainV03_56_B.Throttle_m.Throttle4 = mainV03_56_B.Throttle4_j;
    mainV03_56_B.Throttle_m.Throttle5 = mainV03_56_B.Throttle5_a;
  }

  /* End of Outputs for SubSystem: '<S3>/Fixed-Wing - Cruise' */

  /* Product: '<S3>/Product1' */
  mainV03_56_B.Product1_m0[0] = (real_T)rtb_Compare_a *
    mainV03_56_B.Actuators.deltae;
  mainV03_56_B.Product1_m0[1] = (real_T)rtb_Compare_a *
    mainV03_56_B.Actuators.deltar;
  mainV03_56_B.Product1_m0[2] = (real_T)rtb_Compare_a *
    mainV03_56_B.Actuators.deltafr;
  mainV03_56_B.Product1_m0[3] = (real_T)rtb_Compare_a *
    mainV03_56_B.Actuators.deltafl;

  /* RelationalOperator: '<S232>/Compare' incorporates:
   *  Constant: '<S232>/Constant'
   */
  mainV03_56_B.Compare_a = (mainV03_56_B.DataTypeConversion ==
    mainV03_56_P.CompareToConstant4_const);
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* SignalConversion: '<S3>/HiddenBuf_InsertedFor_TakeOff_at_inport_2' */
    mainV03_56_B.HiddenBuf_InsertedFor_TakeOff_at_inport_2 =
      mainV03_56_B.Compare_a;

    /* Outputs for Enabled SubSystem: '<S3>/TakeOff' incorporates:
     *  EnablePort: '<S168>/Enable'
     */
    if (rtmIsMajorTimeStep(mainV03_56_M)) {
      if (mainV03_56_B.HiddenBuf_InsertedFor_TakeOff_at_inport_2) {
        if (!mainV03_56_DW.TakeOff_MODE) {
          mainV03_56_DW.TakeOff_MODE = true;
        }
      } else {
        if (mainV03_56_DW.TakeOff_MODE) {
          mainV03_56_DW.TakeOff_MODE = false;
        }
      }
    }

    /* End of Outputs for SubSystem: '<S3>/TakeOff' */
  }

  /* Outputs for Enabled SubSystem: '<S3>/TakeOff' incorporates:
   *  EnablePort: '<S168>/Enable'
   */
  if (mainV03_56_DW.TakeOff_MODE) {
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* Gain: '<S168>/Gain' incorporates:
       *  Constant: '<S168>/Constant9'
       */
      mainV03_56_B.actuators.deltae = mainV03_56_P.Gain_Gain_iw *
        mainV03_56_P.deltae_degrees;

      /* Gain: '<S168>/Gain1' incorporates:
       *  Constant: '<S168>/Constant10'
       */
      mainV03_56_B.actuators.deltar = mainV03_56_P.Gain1_Gain_kz *
        mainV03_56_P.deltar_degrees;

      /* Gain: '<S168>/Gain2' incorporates:
       *  Constant: '<S168>/Constant1'
       */
      mainV03_56_B.actuators.deltafr = mainV03_56_P.Gain2_Gain_i *
        mainV03_56_P.deltafr_degrees;

      /* Gain: '<S168>/Gain3' incorporates:
       *  Constant: '<S168>/Constant2'
       */
      mainV03_56_B.actuators.deltafl = mainV03_56_P.Gain3_Gain_c *
        mainV03_56_P.deltafl_degrees;
    }

    /* Sum: '<S168>/Sum1' incorporates:
     *  Constant: '<S168>/Constant15'
     */
    mainV03_56_B.Sum1_b = mainV03_56_P.Constant15_Value_j -
      mainV03_56_B.Sensors.LLAMeas.AltitudeMeas_m;

    /* Saturate: '<S168>/Saturation' */
    if (mainV03_56_B.Sum1_b > mainV03_56_P.Saturation_UpperSat_o) {
      mainV03_56_B.Saturation_c = mainV03_56_P.Saturation_UpperSat_o;
    } else if (mainV03_56_B.Sum1_b < mainV03_56_P.Saturation_LowerSat_j) {
      mainV03_56_B.Saturation_c = mainV03_56_P.Saturation_LowerSat_j;
    } else {
      mainV03_56_B.Saturation_c = mainV03_56_B.Sum1_b;
    }

    /* End of Saturate: '<S168>/Saturation' */

    /* Gain: '<S233>/Proportional Gain' */
    mainV03_56_B.ProportionalGain = mainV03_56_P.PHeight_P_d *
      mainV03_56_B.Saturation_c;

    /* Integrator: '<S233>/Integrator' */
    mainV03_56_B.Integrator = mainV03_56_X.Integrator_CSTATE_p;

    /* Gain: '<S233>/Derivative Gain' */
    mainV03_56_B.DerivativeGain = mainV03_56_P.PHeight_D_n *
      mainV03_56_B.Saturation_c;

    /* Integrator: '<S233>/Filter' */
    mainV03_56_B.Filter = mainV03_56_X.Filter_CSTATE;

    /* Sum: '<S233>/SumD' */
    mainV03_56_B.SumD = mainV03_56_B.DerivativeGain - mainV03_56_B.Filter;

    /* Gain: '<S233>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient = mainV03_56_P.PHeight_N_h *
      mainV03_56_B.SumD;

    /* Sum: '<S233>/Sum' */
    mainV03_56_B.Sum_jd = (mainV03_56_B.ProportionalGain +
      mainV03_56_B.Integrator) + mainV03_56_B.FilterCoefficient;

    /* Saturate: '<S233>/Saturate' */
    if (mainV03_56_B.Sum_jd > mainV03_56_P.PHeight_UpperSaturationLimit_g) {
      mainV03_56_B.Saturate = mainV03_56_P.PHeight_UpperSaturationLimit_g;
    } else if (mainV03_56_B.Sum_jd < mainV03_56_P.PHeight_LowerSaturationLimit_g)
    {
      mainV03_56_B.Saturate = mainV03_56_P.PHeight_LowerSaturationLimit_g;
    } else {
      mainV03_56_B.Saturate = mainV03_56_B.Sum_jd;
    }

    /* End of Saturate: '<S233>/Saturate' */

    /* Sum: '<S168>/Sum12' */
    mainV03_56_B.Sum12 = mainV03_56_B.Saturate - mainV03_56_B.Sensors.alphaMeas;

    /* Gain: '<S235>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_f = mainV03_56_P.PIDController1_P_a *
      mainV03_56_B.Sum12;

    /* Integrator: '<S235>/Integrator' */
    mainV03_56_B.Integrator_l = mainV03_56_X.Integrator_CSTATE_i;

    /* Gain: '<S235>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_d = mainV03_56_P.PIDController1_D_a *
      mainV03_56_B.Sum12;

    /* Integrator: '<S235>/Filter' */
    mainV03_56_B.Filter_h = mainV03_56_X.Filter_CSTATE_e;

    /* Sum: '<S235>/SumD' */
    mainV03_56_B.SumD_f = mainV03_56_B.DerivativeGain_d - mainV03_56_B.Filter_h;

    /* Gain: '<S235>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_l = mainV03_56_P.PIDController1_N_m *
      mainV03_56_B.SumD_f;

    /* Sum: '<S235>/Sum' */
    mainV03_56_B.Sum_pn = (mainV03_56_B.ProportionalGain_f +
      mainV03_56_B.Integrator_l) + mainV03_56_B.FilterCoefficient_l;

    /* Sum: '<S168>/Sum13' */
    mainV03_56_B.PitchRateIn = mainV03_56_B.Sum_pn -
      mainV03_56_B.Sensors.OmegaMeas_body.qMeas;

    /* Gain: '<S234>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_n = mainV03_56_P.PIDController_P_p *
      mainV03_56_B.PitchRateIn;

    /* Integrator: '<S234>/Integrator' */
    mainV03_56_B.Integrator_e = mainV03_56_X.Integrator_CSTATE_n;

    /* Gain: '<S234>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_f = mainV03_56_P.PIDController_D_a *
      mainV03_56_B.PitchRateIn;

    /* Integrator: '<S234>/Filter' */
    mainV03_56_B.Filter_f = mainV03_56_X.Filter_CSTATE_b;

    /* Sum: '<S234>/SumD' */
    mainV03_56_B.SumD_i = mainV03_56_B.DerivativeGain_f - mainV03_56_B.Filter_f;

    /* Gain: '<S234>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_c = mainV03_56_P.PIDController_N_a *
      mainV03_56_B.SumD_i;

    /* Sum: '<S234>/Sum' */
    mainV03_56_B.Sum_g1 = (mainV03_56_B.ProportionalGain_n +
      mainV03_56_B.Integrator_e) + mainV03_56_B.FilterCoefficient_c;

    /* Saturate: '<S234>/Saturate' */
    if (mainV03_56_B.Sum_g1 > mainV03_56_P.PIDController_UpperSaturationLimit_b)
    {
      mainV03_56_B.Saturate_i =
        mainV03_56_P.PIDController_UpperSaturationLimit_b;
    } else if (mainV03_56_B.Sum_g1 <
               mainV03_56_P.PIDController_LowerSaturationLimit_l) {
      mainV03_56_B.Saturate_i =
        mainV03_56_P.PIDController_LowerSaturationLimit_l;
    } else {
      mainV03_56_B.Saturate_i = mainV03_56_B.Sum_g1;
    }

    /* End of Saturate: '<S234>/Saturate' */
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* Gain: '<S168>/Gain4' incorporates:
       *  Constant: '<S168>/Constant'
       */
      mainV03_56_B.deltae_k = mainV03_56_P.Gain4_Gain_fu * mainV03_56_P.delta_e;
    }

    /* ManualSwitch: '<S168>/Manual Switch2' */
    if (mainV03_56_P.ManualSwitch2_CurrentSetting_i == 1) {
      mainV03_56_B.BusCreator1_l.deltae = mainV03_56_B.Saturate_i;
    } else {
      mainV03_56_B.BusCreator1_l.deltae = mainV03_56_B.deltae_k;
    }

    /* End of ManualSwitch: '<S168>/Manual Switch2' */

    /* Sum: '<S168>/Sum3' incorporates:
     *  Constant: '<S168>/Constant14'
     */
    mainV03_56_B.Sum3_j = mainV03_56_P.Constant14_Value_c -
      mainV03_56_B.Sensors.EulerMeas[2];

    /* Gain: '<S237>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_l = mainV03_56_P.PIDController3_P_a *
      mainV03_56_B.Sum3_j;

    /* Integrator: '<S237>/Integrator' */
    mainV03_56_B.Integrator_f = mainV03_56_X.Integrator_CSTATE_d;

    /* Gain: '<S237>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_a = mainV03_56_P.PIDController3_D_g *
      mainV03_56_B.Sum3_j;

    /* Integrator: '<S237>/Filter' */
    mainV03_56_B.Filter_d = mainV03_56_X.Filter_CSTATE_g;

    /* Sum: '<S237>/SumD' */
    mainV03_56_B.SumD_ig = mainV03_56_B.DerivativeGain_a - mainV03_56_B.Filter_d;

    /* Gain: '<S237>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_k = mainV03_56_P.PIDController3_N_h *
      mainV03_56_B.SumD_ig;

    /* Sum: '<S237>/Sum' */
    mainV03_56_B.Sum_jk = (mainV03_56_B.ProportionalGain_l +
      mainV03_56_B.Integrator_f) + mainV03_56_B.FilterCoefficient_k;

    /* Sum: '<S168>/Sum2' */
    mainV03_56_B.Sum2_i = mainV03_56_B.Sum_jk -
      mainV03_56_B.Sensors.OmegaMeas_body.rMeas;

    /* Gain: '<S236>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_g = mainV03_56_P.PIDController2_P_p *
      mainV03_56_B.Sum2_i;

    /* Integrator: '<S236>/Integrator' */
    mainV03_56_B.Integrator_g = mainV03_56_X.Integrator_CSTATE_k;

    /* Gain: '<S236>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_h = mainV03_56_P.PIDController2_D_b *
      mainV03_56_B.Sum2_i;

    /* Integrator: '<S236>/Filter' */
    mainV03_56_B.Filter_dv = mainV03_56_X.Filter_CSTATE_gp;

    /* Sum: '<S236>/SumD' */
    mainV03_56_B.SumD_c = mainV03_56_B.DerivativeGain_h - mainV03_56_B.Filter_dv;

    /* Gain: '<S236>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_o = mainV03_56_P.PIDController2_N_m *
      mainV03_56_B.SumD_c;
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* Gain: '<S168>/Gain5' incorporates:
       *  Constant: '<S168>/Constant11'
       */
      mainV03_56_B.deltar_m = mainV03_56_P.Gain5_Gain_k *
        mainV03_56_P.Constant11_Value_d;
    }

    /* ManualSwitch: '<S168>/Manual Switch3' */
    if (mainV03_56_P.ManualSwitch3_CurrentSetting_p == 1) {
      /* Sum: '<S236>/Sum' */
      mainV03_56_B.Sum_jn = (mainV03_56_B.ProportionalGain_g +
        mainV03_56_B.Integrator_g) + mainV03_56_B.FilterCoefficient_o;
      mainV03_56_B.BusCreator1_l.deltar = mainV03_56_B.Sum_jn;
    } else {
      mainV03_56_B.BusCreator1_l.deltar = mainV03_56_B.deltar_m;
    }

    /* End of ManualSwitch: '<S168>/Manual Switch3' */

    /* Sum: '<S168>/Sum5' incorporates:
     *  Constant: '<S168>/Constant16'
     */
    mainV03_56_B.Sum5_h = mainV03_56_P.Constant16_Value_m -
      mainV03_56_B.Sensors.EulerMeas[0];

    /* Gain: '<S239>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_o = mainV03_56_P.PIDController5_P_b *
      mainV03_56_B.Sum5_h;

    /* Integrator: '<S239>/Integrator' */
    mainV03_56_B.Integrator_lp = mainV03_56_X.Integrator_CSTATE_m;

    /* Gain: '<S239>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_e = mainV03_56_P.PIDController5_D_e *
      mainV03_56_B.Sum5_h;

    /* Integrator: '<S239>/Filter' */
    mainV03_56_B.Filter_fj = mainV03_56_X.Filter_CSTATE_k;

    /* Sum: '<S239>/SumD' */
    mainV03_56_B.SumD_l = mainV03_56_B.DerivativeGain_e - mainV03_56_B.Filter_fj;

    /* Gain: '<S239>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_km = mainV03_56_P.PIDController5_N_m *
      mainV03_56_B.SumD_l;

    /* Sum: '<S239>/Sum' */
    mainV03_56_B.Sum_ck = (mainV03_56_B.ProportionalGain_o +
      mainV03_56_B.Integrator_lp) + mainV03_56_B.FilterCoefficient_km;

    /* Sum: '<S168>/Sum4' */
    mainV03_56_B.Sum4_h = mainV03_56_B.Sum_ck -
      mainV03_56_B.Sensors.OmegaMeas_body.pMeas;

    /* Gain: '<S238>/Proportional Gain' */
    mainV03_56_B.ProportionalGain_a = mainV03_56_P.PIDController4_P_o *
      mainV03_56_B.Sum4_h;

    /* Integrator: '<S238>/Integrator' */
    mainV03_56_B.Integrator_lx = mainV03_56_X.Integrator_CSTATE_py;

    /* Gain: '<S238>/Derivative Gain' */
    mainV03_56_B.DerivativeGain_j = mainV03_56_P.PIDController4_D_nv *
      mainV03_56_B.Sum4_h;

    /* Integrator: '<S238>/Filter' */
    mainV03_56_B.Filter_dd = mainV03_56_X.Filter_CSTATE_n;

    /* Sum: '<S238>/SumD' */
    mainV03_56_B.SumD_h = mainV03_56_B.DerivativeGain_j - mainV03_56_B.Filter_dd;

    /* Gain: '<S238>/Filter Coefficient' */
    mainV03_56_B.FilterCoefficient_m = mainV03_56_P.PIDController4_N_f *
      mainV03_56_B.SumD_h;

    /* Sum: '<S238>/Sum' */
    mainV03_56_B.Sum_f0 = (mainV03_56_B.ProportionalGain_a +
      mainV03_56_B.Integrator_lx) + mainV03_56_B.FilterCoefficient_m;
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* Gain: '<S168>/Gain6' incorporates:
       *  Constant: '<S168>/Constant12'
       */
      mainV03_56_B.deltafr_j = mainV03_56_P.Gain6_Gain_d *
        mainV03_56_P.Constant12_Value_g;

      /* Gain: '<S168>/Gain7' incorporates:
       *  Constant: '<S168>/Constant13'
       */
      mainV03_56_B.deltafl_o = mainV03_56_P.Gain7_Gain_e *
        mainV03_56_P.Constant13_Value_c;
    }

    /* ManualSwitch: '<S168>/Manual Switch4' */
    if (mainV03_56_P.ManualSwitch4_CurrentSetting_e == 1) {
      /* Sum: '<S168>/Sum6' incorporates:
       *  Constant: '<S168>/Constant17'
       */
      mainV03_56_B.Sum6 = mainV03_56_B.Sum_f0 + mainV03_56_P.Constant17_Value_b;
      mainV03_56_B.BusCreator1_l.deltafr = mainV03_56_B.Sum6;
    } else {
      mainV03_56_B.BusCreator1_l.deltafr = mainV03_56_B.deltafr_j;
    }

    /* End of ManualSwitch: '<S168>/Manual Switch4' */

    /* ManualSwitch: '<S168>/Manual Switch6' */
    if (mainV03_56_P.ManualSwitch6_CurrentSetting_k == 1) {
      /* Gain: '<S168>/Gain8' */
      mainV03_56_B.Gain8 = mainV03_56_P.Gain8_Gain_o * mainV03_56_B.Sum_f0;

      /* Sum: '<S168>/Sum' incorporates:
       *  Constant: '<S168>/Constant17'
       */
      mainV03_56_B.Sum_hw = mainV03_56_B.Gain8 + mainV03_56_P.Constant17_Value_b;
      mainV03_56_B.BusCreator1_l.deltafl = mainV03_56_B.Sum_hw;
    } else {
      mainV03_56_B.BusCreator1_l.deltafl = mainV03_56_B.deltafl_o;
    }

    /* End of ManualSwitch: '<S168>/Manual Switch6' */

    /* ManualSwitch: '<S168>/Manual Switch1' */
    if (mainV03_56_P.ManualSwitch1_CurrentSetting_o == 1) {
      mainV03_56_B.ManualSwitch1 = mainV03_56_B.BusCreator1_l;
    } else {
      mainV03_56_B.ManualSwitch1 = mainV03_56_B.actuators;
    }

    /* End of ManualSwitch: '<S168>/Manual Switch1' */

    /* Gain: '<S233>/Integral Gain' */
    mainV03_56_B.IntegralGain = mainV03_56_P.PHeight_I_f *
      mainV03_56_B.Saturation_c;

    /* Gain: '<S234>/Integral Gain' */
    mainV03_56_B.IntegralGain_d = mainV03_56_P.PIDController_I_c *
      mainV03_56_B.PitchRateIn;

    /* Gain: '<S235>/Integral Gain' */
    mainV03_56_B.IntegralGain_i = mainV03_56_P.PIDController1_I_n *
      mainV03_56_B.Sum12;

    /* Gain: '<S236>/Integral Gain' */
    mainV03_56_B.IntegralGain_f = mainV03_56_P.PIDController2_I_g *
      mainV03_56_B.Sum2_i;

    /* Gain: '<S237>/Integral Gain' */
    mainV03_56_B.IntegralGain_a = mainV03_56_P.PIDController3_I_f *
      mainV03_56_B.Sum3_j;

    /* Gain: '<S238>/Integral Gain' */
    mainV03_56_B.IntegralGain_j = mainV03_56_P.PIDController4_I_o *
      mainV03_56_B.Sum4_h;

    /* Gain: '<S239>/Integral Gain' */
    mainV03_56_B.IntegralGain_b = mainV03_56_P.PIDController5_I_b *
      mainV03_56_B.Sum5_h;
    if (rtmIsMajorTimeStep(mainV03_56_M) &&
        mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
      /* ManualSwitch: '<S168>/Manual Switch5' incorporates:
       *  BusCreator: '<S168>/Bus Creator3'
       *  Constant: '<S168>/Constant3'
       *  Constant: '<S168>/Constant4'
       *  Constant: '<S168>/Constant5'
       *  Constant: '<S168>/Constant7'
       *  Constant: '<S168>/Constant8'
       *  Constant: '<S168>/deltat_10'
       *  Constant: '<S168>/deltat_6'
       *  Constant: '<S168>/deltat_7'
       *  Constant: '<S168>/deltat_8'
       *  Constant: '<S168>/deltat_9'
       */
      if (mainV03_56_P.ManualSwitch5_CurrentSetting_h == 1) {
        mainV03_56_B.ManualSwitch5.Throttle1 = mainV03_56_P.Constant4_Value_g;
        mainV03_56_B.ManualSwitch5.Throttle2 = mainV03_56_P.Constant5_Value_p;
        mainV03_56_B.ManualSwitch5.Throttle3 = mainV03_56_P.Constant7_Value_c;
        mainV03_56_B.ManualSwitch5.Throttle4 = mainV03_56_P.Constant8_Value_o;
        mainV03_56_B.ManualSwitch5.Throttle5 = mainV03_56_P.Constant3_Value_l;
      } else {
        mainV03_56_B.ManualSwitch5.Throttle1 = mainV03_56_P.deltat_1;
        mainV03_56_B.ManualSwitch5.Throttle2 = mainV03_56_P.deltat_2;
        mainV03_56_B.ManualSwitch5.Throttle3 = mainV03_56_P.deltat_3;
        mainV03_56_B.ManualSwitch5.Throttle4 = mainV03_56_P.deltat_4;
        mainV03_56_B.ManualSwitch5.Throttle5 = mainV03_56_P.deltat_5;
      }

      /* End of ManualSwitch: '<S168>/Manual Switch5' */
    }
  }

  /* End of Outputs for SubSystem: '<S3>/TakeOff' */

  /* Product: '<S3>/Product9' */
  mainV03_56_B.Product9[0] = (real_T)mainV03_56_B.Compare_a *
    mainV03_56_B.ManualSwitch1.deltae;
  mainV03_56_B.Product9[1] = (real_T)mainV03_56_B.Compare_a *
    mainV03_56_B.ManualSwitch1.deltar;
  mainV03_56_B.Product9[2] = (real_T)mainV03_56_B.Compare_a *
    mainV03_56_B.ManualSwitch1.deltafr;
  mainV03_56_B.Product9[3] = (real_T)mainV03_56_B.Compare_a *
    mainV03_56_B.ManualSwitch1.deltafl;

  /* Sum: '<S3>/Sum1' */
  mainV03_56_B.ActuatorsCmd_b[0] = (((mainV03_56_B.Product7[0] +
    mainV03_56_B.Product5[0]) + mainV03_56_B.Product3[0]) +
    mainV03_56_B.Product1_m0[0]) + mainV03_56_B.Product9[0];
  mainV03_56_B.ActuatorsCmd_b[1] = (((mainV03_56_B.Product7[1] +
    mainV03_56_B.Product5[1]) + mainV03_56_B.Product3[1]) +
    mainV03_56_B.Product1_m0[1]) + mainV03_56_B.Product9[1];
  mainV03_56_B.ActuatorsCmd_b[2] = (((mainV03_56_B.Product7[2] +
    mainV03_56_B.Product5[2]) + mainV03_56_B.Product3[2]) +
    mainV03_56_B.Product1_m0[2]) + mainV03_56_B.Product9[2];
  mainV03_56_B.ActuatorsCmd_b[3] = (((mainV03_56_B.Product7[3] +
    mainV03_56_B.Product5[3]) + mainV03_56_B.Product3[3]) +
    mainV03_56_B.Product1_m0[3]) + mainV03_56_B.Product9[3];

  /* BusCreator: '<S3>/Bus Creator' */
  mainV03_56_B.ActuatorsCmd.deltae = mainV03_56_B.ActuatorsCmd_b[0];
  mainV03_56_B.ActuatorsCmd.deltar = mainV03_56_B.ActuatorsCmd_b[1];
  mainV03_56_B.ActuatorsCmd.deltafr = mainV03_56_B.ActuatorsCmd_b[2];
  mainV03_56_B.ActuatorsCmd.deltafl = mainV03_56_B.ActuatorsCmd_b[3];

  /* Product: '<S3>/Product6' */
  mainV03_56_B.Product6[0] = (real_T)mainV03_56_B.Compare *
    mainV03_56_B.ManualSwitch3_i.Throttle1;
  mainV03_56_B.Product6[1] = (real_T)mainV03_56_B.Compare *
    mainV03_56_B.ManualSwitch3_i.Throttle2;
  mainV03_56_B.Product6[2] = (real_T)mainV03_56_B.Compare *
    mainV03_56_B.ManualSwitch3_i.Throttle3;
  mainV03_56_B.Product6[3] = (real_T)mainV03_56_B.Compare *
    mainV03_56_B.ManualSwitch3_i.Throttle4;
  mainV03_56_B.Product6[4] = (real_T)mainV03_56_B.Compare *
    mainV03_56_B.ManualSwitch3_i.Throttle5;

  /* Product: '<S3>/Product4' */
  mainV03_56_B.Product4_o[0] = (real_T)mainV03_56_B.Compare_j *
    mainV03_56_B.ManualSwitch3.Throttle1;
  mainV03_56_B.Product4_o[1] = (real_T)mainV03_56_B.Compare_j *
    mainV03_56_B.ManualSwitch3.Throttle2;
  mainV03_56_B.Product4_o[2] = (real_T)mainV03_56_B.Compare_j *
    mainV03_56_B.ManualSwitch3.Throttle3;
  mainV03_56_B.Product4_o[3] = (real_T)mainV03_56_B.Compare_j *
    mainV03_56_B.ManualSwitch3.Throttle4;
  mainV03_56_B.Product4_o[4] = (real_T)mainV03_56_B.Compare_j *
    mainV03_56_B.ManualSwitch3.Throttle5;

  /* Product: '<S3>/Product2' */
  mainV03_56_B.Product2_a[0] = (real_T)mainV03_56_B.Compare_g *
    mainV03_56_B.ManualSwitch5_j.Throttle1;
  mainV03_56_B.Product2_a[1] = (real_T)mainV03_56_B.Compare_g *
    mainV03_56_B.ManualSwitch5_j.Throttle2;
  mainV03_56_B.Product2_a[2] = (real_T)mainV03_56_B.Compare_g *
    mainV03_56_B.ManualSwitch5_j.Throttle3;
  mainV03_56_B.Product2_a[3] = (real_T)mainV03_56_B.Compare_g *
    mainV03_56_B.ManualSwitch5_j.Throttle4;
  mainV03_56_B.Product2_a[4] = (real_T)mainV03_56_B.Compare_g *
    mainV03_56_B.ManualSwitch5_j.Throttle5;

  /* Product: '<S3>/Product' */
  mainV03_56_B.Product_ki[0] = (real_T)rtb_Compare_a *
    mainV03_56_B.Throttle_m.Throttle1;
  mainV03_56_B.Product_ki[1] = (real_T)rtb_Compare_a *
    mainV03_56_B.Throttle_m.Throttle2;
  mainV03_56_B.Product_ki[2] = (real_T)rtb_Compare_a *
    mainV03_56_B.Throttle_m.Throttle3;
  mainV03_56_B.Product_ki[3] = (real_T)rtb_Compare_a *
    mainV03_56_B.Throttle_m.Throttle4;
  mainV03_56_B.Product_ki[4] = (real_T)rtb_Compare_a *
    mainV03_56_B.Throttle_m.Throttle5;

  /* Product: '<S3>/Product8' */
  mainV03_56_B.Product8[0] = (real_T)mainV03_56_B.Compare_a *
    mainV03_56_B.ManualSwitch5.Throttle1;
  mainV03_56_B.Product8[1] = (real_T)mainV03_56_B.Compare_a *
    mainV03_56_B.ManualSwitch5.Throttle2;
  mainV03_56_B.Product8[2] = (real_T)mainV03_56_B.Compare_a *
    mainV03_56_B.ManualSwitch5.Throttle3;
  mainV03_56_B.Product8[3] = (real_T)mainV03_56_B.Compare_a *
    mainV03_56_B.ManualSwitch5.Throttle4;
  mainV03_56_B.Product8[4] = (real_T)mainV03_56_B.Compare_a *
    mainV03_56_B.ManualSwitch5.Throttle5;

  /* Sum: '<S3>/Sum2' */
  for (s155_iter = 0; s155_iter < 5; s155_iter++) {
    mainV03_56_B.Throttle_i[s155_iter] = (((mainV03_56_B.Product6[s155_iter] +
      mainV03_56_B.Product4_o[s155_iter]) + mainV03_56_B.Product2_a[s155_iter])
      + mainV03_56_B.Product_ki[s155_iter]) + mainV03_56_B.Product8[s155_iter];
  }

  /* End of Sum: '<S3>/Sum2' */

  /* BusCreator: '<S3>/Bus Creator1' */
  mainV03_56_B.Throttle.Throttle1 = mainV03_56_B.Throttle_i[0];
  mainV03_56_B.Throttle.Throttle2 = mainV03_56_B.Throttle_i[1];
  mainV03_56_B.Throttle.Throttle3 = mainV03_56_B.Throttle_i[2];
  mainV03_56_B.Throttle.Throttle4 = mainV03_56_B.Throttle_i[3];
  mainV03_56_B.Throttle.Throttle5 = mainV03_56_B.Throttle_i[4];

  /* Sum: '<S242>/Sum' */
  mainV03_56_B.Sum_ba[0] = mainV03_56_B.plantData.V_body[0] -
    mainV03_56_B.envData.windVelocity[0];
  mainV03_56_B.Sum_ba[1] = mainV03_56_B.plantData.V_body[1] -
    mainV03_56_B.envData.windVelocity[1];
  mainV03_56_B.Sum_ba[2] = mainV03_56_B.plantData.V_body[2] -
    mainV03_56_B.envData.windVelocity[2];

  /* Trigonometry: '<S261>/Incidence' */
  rtb_Rn = rt_atan2d_snf(mainV03_56_B.Sum_ba[2], mainV03_56_B.Sum_ba[0]);

  /* Trigonometry: '<S259>/Trigonometric Function2' */
  rtb_UnaryMinus_j = sin(rtb_Rn);
  rtb_UnaryMinus_o = cos(rtb_Rn);

  /* SignalConversion: '<S286>/ConcatBufferAtVector ConcatenateIn1' */
  mainV03_56_B.VectorConcatenate_l[0] = rtb_UnaryMinus_o;
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Constant: '<S259>/Constant1' */
    mainV03_56_B.VectorConcatenate_l[1] = mainV03_56_P.Constant1_Value_i;
  }

  /* UnaryMinus: '<S259>/Unary Minus' */
  mainV03_56_B.VectorConcatenate_l[2] = -rtb_UnaryMinus_j;
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Constant: '<S259>/Constant' */
    mainV03_56_B.VectorConcatenate_l[3] = mainV03_56_P.Constant_Value_bw;

    /* Constant: '<S259>/Constant2' */
    mainV03_56_B.VectorConcatenate_l[4] = mainV03_56_P.Constant2_Value_ne;

    /* Constant: '<S259>/Constant4' */
    mainV03_56_B.VectorConcatenate_l[5] = mainV03_56_P.Constant4_Value_p1;
  }

  /* SignalConversion: '<S286>/ConcatBufferAtVector ConcatenateIn7' */
  mainV03_56_B.VectorConcatenate_l[6] = rtb_UnaryMinus_j;
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Constant: '<S259>/Constant3' */
    mainV03_56_B.VectorConcatenate_l[7] = mainV03_56_P.Constant3_Value_p;
  }

  /* SignalConversion: '<S286>/ConcatBufferAtVector ConcatenateIn9' */
  mainV03_56_B.VectorConcatenate_l[8] = rtb_UnaryMinus_o;

  /* PreLookup: '<S242>/preAlpha1' */
  mainV03_56_B.alphaBus_g.idxalpha = plook_s32dd_bincpa
    (mainV03_56_B.plantData.alpha, mainV03_56_P.preAlpha1_BreakpointsData, 8U,
     &mainV03_56_B.alphaBus_g.falpha, &mainV03_56_DW.preAlpha1_DWORK1);
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* BusCreator: '<S242>/Bus Creator17' incorporates:
     *  Constant: '<S242>/Constant '
     *  Constant: '<S242>/Constant 1'
     */
    mainV03_56_B.alphaBus_a.idxalpha = mainV03_56_P.Constant_Value_g0l;
    mainV03_56_B.alphaBus_a.falpha = mainV03_56_P.Constant1_Value_dt;
  }

  /* Switch: '<S242>/Switch7' incorporates:
   *  Constant: '<S242>/LD.flags.alpha.length'
   */
  if (mainV03_56_P.LDflagsalphalength_Value > mainV03_56_P.Switch7_Threshold) {
    mainV03_56_B.analysisCasesBus_j.alpha = mainV03_56_B.alphaBus_g;
  } else {
    mainV03_56_B.analysisCasesBus_j.alpha = mainV03_56_B.alphaBus_a;
  }

  /* End of Switch: '<S242>/Switch7' */

  /* PreLookup: '<S242>/preBeta' */
  mainV03_56_B.betaBus_p.idxbeta = plook_s32dd_bincpa
    (mainV03_56_B.plantData.beta, mainV03_56_P.preBeta_BreakpointsData, 1U,
     &mainV03_56_B.betaBus_p.fbeta, &mainV03_56_DW.preBeta_DWORK1);
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* BusCreator: '<S242>/Bus Creator2' incorporates:
     *  Constant: '<S242>/Constant 14'
     *  Constant: '<S242>/Constant 15'
     */
    mainV03_56_B.betaBus_g.idxbeta = mainV03_56_P.Constant14_Value_b;
    mainV03_56_B.betaBus_g.fbeta = mainV03_56_P.Constant15_Value_js;
  }

  /* Switch: '<S242>/Switch1' incorporates:
   *  Constant: '<S242>/LD.flags.beta.length'
   */
  if (mainV03_56_P.LDflagsbetalength_Value > mainV03_56_P.Switch1_Threshold_f) {
    mainV03_56_B.analysisCasesBus_j.beta = mainV03_56_B.betaBus_p;
  } else {
    mainV03_56_B.analysisCasesBus_j.beta = mainV03_56_B.betaBus_g;
  }

  /* End of Switch: '<S242>/Switch1' */

  /* PreLookup: '<S242>/h  look-up1' */
  mainV03_56_B.altBus_n.idxalt = plook_s32dd_bincpa
    (mainV03_56_B.plantData.LLA.Altitude_m,
     mainV03_56_P.hlookup1_BreakpointsData, 1U, &mainV03_56_B.altBus_n.falt,
     &mainV03_56_DW.hlookup1_DWORK1);
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* BusCreator: '<S242>/Bus Creator18' incorporates:
     *  Constant: '<S242>/Constant 2'
     *  Constant: '<S242>/Constant 3'
     */
    mainV03_56_B.altBus_b.idxalt = mainV03_56_P.Constant2_Value_am;
    mainV03_56_B.altBus_b.falt = mainV03_56_P.Constant3_Value_b;
  }

  /* Switch: '<S242>/Switch8' incorporates:
   *  Constant: '<S242>/LD.flags.alt.length'
   */
  if (mainV03_56_P.LDflagsaltlength_Value > mainV03_56_P.Switch8_Threshold) {
    mainV03_56_B.analysisCasesBus_j.alt = mainV03_56_B.altBus_n;
  } else {
    mainV03_56_B.analysisCasesBus_j.alt = mainV03_56_B.altBus_b;
  }

  /* End of Switch: '<S242>/Switch8' */

  /* PreLookup: '<S242>/ look-up1' */
  mainV03_56_B.xcgBus_o.idxxcg = plook_s32dd_bincpa(mainV03_56_B.plantData.CG[0],
    mainV03_56_P.lookup1_BreakpointsData, 1U, &mainV03_56_B.xcgBus_o.fxcg,
    &mainV03_56_DW.lookup1_DWORK1);
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* BusCreator: '<S242>/Bus Creator23' incorporates:
     *  Constant: '<S242>/Constant 4'
     *  Constant: '<S242>/Constant 5'
     */
    mainV03_56_B.xcgBus_p.idxxcg = mainV03_56_P.Constant4_Value_m;
    mainV03_56_B.xcgBus_p.fxcg = mainV03_56_P.Constant5_Value_n;
  }

  /* Switch: '<S242>/Switch9' incorporates:
   *  Constant: '<S242>/LD.flags.xcg.length'
   */
  if (mainV03_56_P.LDflagsxcglength_Value > mainV03_56_P.Switch9_Threshold) {
    mainV03_56_B.analysisCasesBus_j.xcg = mainV03_56_B.xcgBus_o;
  } else {
    mainV03_56_B.analysisCasesBus_j.xcg = mainV03_56_B.xcgBus_p;
  }

  /* End of Switch: '<S242>/Switch9' */

  /* PreLookup: '<S242>/p1' */
  mainV03_56_B.deltaeBus_h.idxdeltae = plook_s32dd_bincpa
    (mainV03_56_B.ActuatorsCmd.deltae, mainV03_56_P.p1_BreakpointsData, 1U,
     &mainV03_56_B.deltaeBus_h.fdeltae, &mainV03_56_DW.p1_DWORK1);
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* BusCreator: '<S242>/Bus Creator19' incorporates:
     *  Constant: '<S242>/Constant 6'
     *  Constant: '<S242>/Constant 7'
     */
    mainV03_56_B.deltaeBus_f.idxdeltae = mainV03_56_P.Constant6_Value_jg;
    mainV03_56_B.deltaeBus_f.fdeltae = mainV03_56_P.Constant7_Value_o;
  }

  /* Switch: '<S242>/Switch10' incorporates:
   *  Constant: '<S242>/LD.flags.deltae.length'
   */
  if (mainV03_56_P.LDflagsdeltaelength_Value > mainV03_56_P.Switch10_Threshold)
  {
    mainV03_56_B.analysisCasesBus_j.deltae = mainV03_56_B.deltaeBus_h;
  } else {
    mainV03_56_B.analysisCasesBus_j.deltae = mainV03_56_B.deltaeBus_f;
  }

  /* End of Switch: '<S242>/Switch10' */

  /* PreLookup: '<S242>/look-up1' */
  mainV03_56_B.deltarBus_a.idxdeltar = plook_s32dd_bincpa
    (mainV03_56_B.ActuatorsCmd.deltar, mainV03_56_P.lookup1_BreakpointsData_g,
     1U, &mainV03_56_B.deltarBus_a.fdeltar, &mainV03_56_DW.lookup1_DWORK1_l);
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* BusCreator: '<S242>/Bus Creator20' incorporates:
     *  Constant: '<S242>/Constant 8'
     *  Constant: '<S242>/Constant 9'
     */
    mainV03_56_B.deltarBus_b.idxdeltar = mainV03_56_P.Constant8_Value_a;
    mainV03_56_B.deltarBus_b.fdeltar = mainV03_56_P.Constant9_Value_a;
  }

  /* Switch: '<S242>/Switch11' incorporates:
   *  Constant: '<S242>/LD.flags.deltar.length'
   */
  if (mainV03_56_P.LDflagsdeltarlength_Value > mainV03_56_P.Switch11_Threshold)
  {
    mainV03_56_B.analysisCasesBus_j.deltar = mainV03_56_B.deltarBus_a;
  } else {
    mainV03_56_B.analysisCasesBus_j.deltar = mainV03_56_B.deltarBus_b;
  }

  /* End of Switch: '<S242>/Switch11' */

  /* PreLookup: '<S242>/fr look-up1' */
  mainV03_56_B.deltafrBus_d.idxdeltafr = plook_s32dd_bincpa
    (mainV03_56_B.ActuatorsCmd.deltafr, mainV03_56_P.frlookup1_BreakpointsData,
     1U, &mainV03_56_B.deltafrBus_d.fdeltafr, &mainV03_56_DW.frlookup1_DWORK1);
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* BusCreator: '<S242>/Bus Creator21' incorporates:
     *  Constant: '<S242>/Constant 10'
     *  Constant: '<S242>/Constant 11'
     */
    mainV03_56_B.deltafrBus_l.idxdeltafr = mainV03_56_P.Constant10_Value_d;
    mainV03_56_B.deltafrBus_l.fdeltafr = mainV03_56_P.Constant11_Value_dd;
  }

  /* Switch: '<S242>/Switch12' incorporates:
   *  Constant: '<S242>/LD.flags.deltafr.length'
   */
  if (mainV03_56_P.LDflagsdeltafrlength_Value > mainV03_56_P.Switch12_Threshold)
  {
    mainV03_56_B.analysisCasesBus_j.deltafr = mainV03_56_B.deltafrBus_d;
  } else {
    mainV03_56_B.analysisCasesBus_j.deltafr = mainV03_56_B.deltafrBus_l;
  }

  /* End of Switch: '<S242>/Switch12' */

  /* PreLookup: '<S242>/pre look-ups1' */
  mainV03_56_B.deltaflBus_i.idxdeltafl = plook_s32dd_bincpa
    (mainV03_56_B.ActuatorsCmd.deltafl, mainV03_56_P.prelookups1_BreakpointsData,
     1U, &mainV03_56_B.deltaflBus_i.fdeltafl, &mainV03_56_DW.prelookups1_DWORK1);
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* BusCreator: '<S242>/Bus Creator22' incorporates:
     *  Constant: '<S242>/Constant 12'
     *  Constant: '<S242>/Constant 13'
     */
    mainV03_56_B.deltaflBus_k.idxdeltafl = mainV03_56_P.Constant12_Value_i;
    mainV03_56_B.deltaflBus_k.fdeltafl = mainV03_56_P.Constant13_Value_cc;
  }

  /* Switch: '<S242>/Switch13' incorporates:
   *  Constant: '<S242>/LD.flags.deltafl.length'
   */
  if (mainV03_56_P.LDflagsdeltafllength_Value > mainV03_56_P.Switch13_Threshold)
  {
    mainV03_56_B.analysisCasesBus_j.deltafl = mainV03_56_B.deltaflBus_i;
  } else {
    mainV03_56_B.analysisCasesBus_j.deltafl = mainV03_56_B.deltaflBus_k;
  }

  /* End of Switch: '<S242>/Switch13' */

  /* Interpolation_n-D: '<S253>/CD0' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_0[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_0[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_0[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_0[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_0[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_0[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_0[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_0[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_0[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_0[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_0[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_0[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_0[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_0[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_0[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_0[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_0[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_0[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_0[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_0[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_0[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_0[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_0[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_0[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_0[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_0[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_0[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_0[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_0[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_0[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_0[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_0[7] = 0;
  }

  mainV03_56_B.CD0 = intrp8d_s32dla_pw(bpIndex_0, frac_0, mainV03_56_P.CD0_Table,
    mainV03_56_P.CD0_dimSize, mainV03_56_P.CD0_maxIndex);

  /* End of Interpolation_n-D: '<S253>/CD0' */

  /* Interpolation_n-D: '<S253>/CDalpha' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_1[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_1[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_1[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_1[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_1[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_1[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_1[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_1[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_1[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_1[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_1[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_1[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_1[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_1[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_1[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_1[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_1[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_1[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_1[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_1[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_1[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_1[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_1[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_1[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_1[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_1[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_1[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_1[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_1[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_1[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_1[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_1[7] = 0;
  }

  mainV03_56_B.CDalpha = intrp8d_s32dla_pw(bpIndex_1, frac_1,
    mainV03_56_P.CDalpha_Table, mainV03_56_P.CDalpha_dimSize,
    mainV03_56_P.CDalpha_maxIndex);

  /* End of Interpolation_n-D: '<S253>/CDalpha' */

  /* Product: '<S253>/Product' */
  mainV03_56_B.Product_d = mainV03_56_B.CDalpha * mainV03_56_B.plantData.alpha;

  /* Interpolation_n-D: '<S253>/CDalpha_dot' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_2[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_2[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_2[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_2[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_2[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_2[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_2[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_2[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_2[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_2[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_2[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_2[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_2[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_2[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_2[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_2[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_2[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_2[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_2[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_2[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_2[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_2[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_2[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_2[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_2[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_2[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_2[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_2[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_2[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_2[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_2[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_2[7] = 0;
  }

  mainV03_56_B.CDalpha_dot = intrp8d_s32dla_pw(bpIndex_2, frac_2,
    mainV03_56_P.CDalpha_dot_Table, mainV03_56_P.CDalpha_dot_dimSize,
    mainV03_56_P.CDalpha_dot_maxIndex);

  /* End of Interpolation_n-D: '<S253>/CDalpha_dot' */

  /* Product: '<S253>/Product1' */
  mainV03_56_B.Product1_o = mainV03_56_B.CDalpha_dot *
    mainV03_56_B.plantData.alpha_dot;

  /* Interpolation_n-D: '<S253>/CDq' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_3[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_3[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_3[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_3[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_3[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_3[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_3[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_3[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_3[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_3[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_3[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_3[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_3[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_3[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_3[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_3[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_3[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_3[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_3[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_3[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_3[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_3[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_3[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_3[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_3[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_3[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_3[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_3[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_3[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_3[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_3[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_3[7] = 0;
  }

  mainV03_56_B.CDq = intrp8d_s32dla_pw(bpIndex_3, frac_3, mainV03_56_P.CDq_Table,
    mainV03_56_P.CDq_dimSize, mainV03_56_P.CDq_maxIndex);

  /* End of Interpolation_n-D: '<S253>/CDq' */

  /* Product: '<S253>/Product2' */
  mainV03_56_B.Product2_m = mainV03_56_B.CDq *
    mainV03_56_B.plantData.Omega_body.q;

  /* Interpolation_n-D: '<S253>/CDdeltae' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_4[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_4[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_4[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_4[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_4[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_4[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_4[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_4[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_4[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_4[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_4[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_4[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_4[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_4[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_4[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_4[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_4[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_4[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_4[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_4[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_4[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_4[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_4[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_4[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_4[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_4[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_4[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_4[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_4[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_4[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_4[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_4[7] = 0;
  }

  mainV03_56_B.CDdeltae = intrp8d_s32dla_pw(bpIndex_4, frac_4,
    mainV03_56_P.CDdeltae_Table, mainV03_56_P.CDdeltae_dimSize,
    mainV03_56_P.CDdeltae_maxIndex);

  /* End of Interpolation_n-D: '<S253>/CDdeltae' */

  /* Product: '<S253>/Product3' */
  mainV03_56_B.Product3_m = mainV03_56_B.CDdeltae *
    mainV03_56_B.ActuatorsCmd.deltae;

  /* Interpolation_n-D: '<S253>/CDdeltafr' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_5[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_5[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_5[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_5[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_5[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_5[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_5[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_5[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_5[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_5[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_5[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_5[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_5[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_5[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_5[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_5[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_5[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_5[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_5[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_5[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_5[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_5[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_5[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_5[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_5[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_5[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_5[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_5[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_5[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_5[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_5[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_5[7] = 0;
  }

  mainV03_56_B.CDdeltafr = intrp8d_s32dla_pw(bpIndex_5, frac_5,
    mainV03_56_P.CDdeltafr_Table, mainV03_56_P.CDdeltafr_dimSize,
    mainV03_56_P.CDdeltafr_maxIndex);

  /* End of Interpolation_n-D: '<S253>/CDdeltafr' */

  /* Product: '<S253>/Product4' */
  mainV03_56_B.Product4_e = mainV03_56_B.CDdeltafr *
    mainV03_56_B.ActuatorsCmd.deltafr;

  /* Interpolation_n-D: '<S253>/CDdeltafrl' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_6[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_6[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_6[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_6[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_6[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_6[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_6[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_6[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_6[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_6[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_6[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_6[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_6[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_6[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_6[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_6[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_6[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_6[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_6[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_6[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_6[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_6[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_6[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_6[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_6[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_6[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_6[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_6[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_6[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_6[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_6[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_6[7] = 0;
  }

  mainV03_56_B.CDdeltafrl = intrp8d_s32dla_pw(bpIndex_6, frac_6,
    mainV03_56_P.CDdeltafrl_Table, mainV03_56_P.CDdeltafrl_dimSize,
    mainV03_56_P.CDdeltafrl_maxIndex);

  /* End of Interpolation_n-D: '<S253>/CDdeltafrl' */

  /* Product: '<S253>/Product5' */
  mainV03_56_B.Product5_a = mainV03_56_B.CDdeltafrl *
    mainV03_56_B.ActuatorsCmd.deltafl;

  /* Sum: '<S253>/Add' */
  mainV03_56_B.Add = (((((mainV03_56_B.CD0 + mainV03_56_B.Product_d) +
    mainV03_56_B.Product1_o) + mainV03_56_B.Product2_m) +
                       mainV03_56_B.Product3_m) + mainV03_56_B.Product4_e) +
    mainV03_56_B.Product5_a;

  /* Interpolation_n-D: '<S255>/CYbeta' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_7[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_7[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_7[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_7[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_7[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_7[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_7[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_7[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_7[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_7[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_7[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_7[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_7[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_7[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_7[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_7[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_7[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_7[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_7[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_7[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_7[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_7[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_7[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_7[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_7[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_7[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_7[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_7[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_7[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_7[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_7[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_7[7] = 0;
  }

  mainV03_56_B.CYbeta = intrp8d_s32dla_pw(bpIndex_7, frac_7,
    mainV03_56_P.CYbeta_Table, mainV03_56_P.CYbeta_dimSize,
    mainV03_56_P.CYbeta_maxIndex);

  /* End of Interpolation_n-D: '<S255>/CYbeta' */

  /* Product: '<S255>/Product' */
  mainV03_56_B.Product_gc = mainV03_56_B.CYbeta * mainV03_56_B.plantData.beta;

  /* Interpolation_n-D: '<S255>/CYp' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_8[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_8[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_8[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_8[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_8[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_8[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_8[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_8[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_8[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_8[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_8[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_8[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_8[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_8[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_8[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_8[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_8[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_8[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_8[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_8[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_8[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_8[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_8[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_8[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_8[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_8[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_8[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_8[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_8[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_8[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_8[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_8[7] = 0;
  }

  mainV03_56_B.CYp = intrp8d_s32dla_pw(bpIndex_8, frac_8, mainV03_56_P.CYp_Table,
    mainV03_56_P.CYp_dimSize, mainV03_56_P.CYp_maxIndex);

  /* End of Interpolation_n-D: '<S255>/CYp' */

  /* Product: '<S255>/Product1' */
  mainV03_56_B.Product1_p = mainV03_56_B.CYp *
    mainV03_56_B.plantData.Omega_body.p;

  /* Interpolation_n-D: '<S255>/CYr' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_9[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_9[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_9[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_9[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_9[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_9[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_9[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_9[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_9[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_9[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_9[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_9[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_9[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_9[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_9[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_9[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_9[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_9[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_9[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_9[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_9[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_9[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_9[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_9[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_9[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_9[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_9[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_9[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_9[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_9[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_9[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_9[7] = 0;
  }

  mainV03_56_B.CYr = intrp8d_s32dla_pw(bpIndex_9, frac_9, mainV03_56_P.CYr_Table,
    mainV03_56_P.CYr_dimSize, mainV03_56_P.CYr_maxIndex);

  /* End of Interpolation_n-D: '<S255>/CYr' */

  /* Product: '<S255>/Product2' */
  mainV03_56_B.Product2_j = mainV03_56_B.CYr *
    mainV03_56_B.plantData.Omega_body.r;

  /* Interpolation_n-D: '<S255>/CYdeltar ' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_a[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_a[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_a[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_a[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_a[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_a[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_a[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_a[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_a[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_a[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_a[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_a[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_a[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_a[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_a[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_a[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_a[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_a[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_a[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_a[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_a[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_a[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_a[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_a[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_a[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_a[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_a[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_a[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_a[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_a[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_a[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_a[7] = 0;
  }

  mainV03_56_B.CYdeltar = intrp8d_s32dla_pw(bpIndex_a, frac_a,
    mainV03_56_P.CYdeltar_Table, mainV03_56_P.CYdeltar_dimSize,
    mainV03_56_P.CYdeltar_maxIndex);

  /* End of Interpolation_n-D: '<S255>/CYdeltar ' */

  /* Product: '<S255>/Product3' */
  mainV03_56_B.Product3_i = mainV03_56_B.CYdeltar *
    mainV03_56_B.ActuatorsCmd.deltar;

  /* Interpolation_n-D: '<S255>/CYdeltafr' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_b[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_b[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_b[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_b[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_b[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_b[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_b[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_b[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_b[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_b[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_b[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_b[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_b[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_b[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_b[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_b[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_b[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_b[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_b[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_b[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_b[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_b[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_b[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_b[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_b[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_b[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_b[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_b[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_b[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_b[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_b[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_b[7] = 0;
  }

  mainV03_56_B.CYdeltafr = intrp8d_s32dla_pw(bpIndex_b, frac_b,
    mainV03_56_P.CYdeltafr_Table, mainV03_56_P.CYdeltafr_dimSize,
    mainV03_56_P.CYdeltafr_maxIndex);

  /* End of Interpolation_n-D: '<S255>/CYdeltafr' */

  /* Product: '<S255>/Product4' */
  mainV03_56_B.Product4_f = mainV03_56_B.CYdeltafr *
    mainV03_56_B.ActuatorsCmd.deltafr;

  /* Interpolation_n-D: '<S255>/CYdeltafl' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_c[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_c[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_c[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_c[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_c[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_c[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_c[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_c[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_c[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_c[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_c[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_c[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_c[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_c[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_c[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_c[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_c[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_c[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_c[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_c[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_c[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_c[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_c[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_c[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_c[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_c[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_c[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_c[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_c[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_c[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_c[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_c[7] = 0;
  }

  mainV03_56_B.CYdeltafl = intrp8d_s32dla_pw(bpIndex_c, frac_c,
    mainV03_56_P.CYdeltafl_Table, mainV03_56_P.CYdeltafl_dimSize,
    mainV03_56_P.CYdeltafl_maxIndex);

  /* End of Interpolation_n-D: '<S255>/CYdeltafl' */

  /* Product: '<S255>/Product5' */
  mainV03_56_B.Product5_f = mainV03_56_B.CYdeltafl *
    mainV03_56_B.ActuatorsCmd.deltafl;

  /* Sum: '<S255>/Add' */
  mainV03_56_B.Add_g = ((((mainV03_56_B.Product_gc + mainV03_56_B.Product1_p) +
    mainV03_56_B.Product2_j) + mainV03_56_B.Product3_i) +
                        mainV03_56_B.Product4_f) + mainV03_56_B.Product5_f;

  /* Interpolation_n-D: '<S254>/CL0' */
  frac_d[0] = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  frac_d[1] = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  frac_d[2] = mainV03_56_B.analysisCasesBus_j.alt.falt;
  frac_d[3] = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  frac_d[4] = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  frac_d[5] = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  frac_d[6] = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  frac_d[7] = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 7) {
    bpIndex_d[0] = 7;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_d[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_d[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 0) {
    bpIndex_d[1] = 0;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_d[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_d[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 0) {
    bpIndex_d[2] = 0;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_d[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_d[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 0) {
    bpIndex_d[3] = 0;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_d[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_d[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 0) {
    bpIndex_d[4] = 0;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_d[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_d[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 0) {
    bpIndex_d[5] = 0;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_d[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_d[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 0) {
    bpIndex_d[6] = 0;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_d[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_d[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 0) {
    bpIndex_d[7] = 0;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_d[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_d[7] = 0;
  }

  mainV03_56_B.CL0 = intrp8d_s32dl_pw(bpIndex_d, frac_d, mainV03_56_P.CL0_Table,
    mainV03_56_P.CL0_dimSize);

  /* End of Interpolation_n-D: '<S254>/CL0' */

  /* Interpolation_n-D: '<S254>/CLalpha' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_e[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_e[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_e[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_e[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_e[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_e[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_e[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_e[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_e[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_e[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_e[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_e[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_e[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_e[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_e[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_e[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_e[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_e[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_e[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_e[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_e[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_e[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_e[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_e[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_e[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_e[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_e[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_e[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_e[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_e[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_e[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_e[7] = 0;
  }

  mainV03_56_B.CLalpha = intrp8d_s32dla_pw(bpIndex_e, frac_e,
    mainV03_56_P.CLalpha_Table, mainV03_56_P.CLalpha_dimSize,
    mainV03_56_P.CLalpha_maxIndex);

  /* End of Interpolation_n-D: '<S254>/CLalpha' */

  /* Product: '<S254>/Product6' */
  mainV03_56_B.Product6_e = mainV03_56_B.CLalpha * mainV03_56_B.plantData.alpha;

  /* Interpolation_n-D: '<S254>/CLalpha_dot' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_f[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_f[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_f[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_f[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_f[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_f[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_f[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_f[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_f[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_f[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_f[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_f[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_f[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_f[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_f[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_f[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_f[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_f[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_f[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_f[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_f[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_f[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_f[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_f[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_f[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_f[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_f[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_f[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_f[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_f[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_f[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_f[7] = 0;
  }

  mainV03_56_B.CLalpha_dot = intrp8d_s32dla_pw(bpIndex_f, frac_f,
    mainV03_56_P.CLalpha_dot_Table, mainV03_56_P.CLalpha_dot_dimSize,
    mainV03_56_P.CLalpha_dot_maxIndex);

  /* End of Interpolation_n-D: '<S254>/CLalpha_dot' */

  /* Product: '<S254>/Product7' */
  mainV03_56_B.Product7_j = mainV03_56_B.CLalpha_dot *
    mainV03_56_B.plantData.alpha_dot;

  /* Interpolation_n-D: '<S254>/CLq' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_g[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_g[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_g[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_g[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_g[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_g[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_g[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_g[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_g[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_g[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_g[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_g[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_g[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_g[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_g[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_g[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_g[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_g[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_g[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_g[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_g[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_g[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_g[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_g[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_g[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_g[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_g[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_g[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_g[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_g[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_g[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_g[7] = 0;
  }

  mainV03_56_B.CLq = intrp8d_s32dla_pw(bpIndex_g, frac_g, mainV03_56_P.CLq_Table,
    mainV03_56_P.CLq_dimSize, mainV03_56_P.CLq_maxIndex);

  /* End of Interpolation_n-D: '<S254>/CLq' */

  /* Product: '<S254>/Product8' */
  mainV03_56_B.Product8_h = mainV03_56_B.CLq *
    mainV03_56_B.plantData.Omega_body.q;

  /* Interpolation_n-D: '<S254>/CLdeltae' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_h[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_h[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_h[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_h[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_h[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_h[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_h[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_h[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_h[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_h[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_h[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_h[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_h[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_h[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_h[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_h[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_h[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_h[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_h[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_h[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_h[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_h[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_h[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_h[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_h[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_h[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_h[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_h[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_h[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_h[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_h[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_h[7] = 0;
  }

  mainV03_56_B.CLdeltae = intrp8d_s32dla_pw(bpIndex_h, frac_h,
    mainV03_56_P.CLdeltae_Table, mainV03_56_P.CLdeltae_dimSize,
    mainV03_56_P.CLdeltae_maxIndex);

  /* End of Interpolation_n-D: '<S254>/CLdeltae' */

  /* Product: '<S254>/Product9' */
  mainV03_56_B.Product9_j = mainV03_56_B.CLdeltae *
    mainV03_56_B.ActuatorsCmd.deltae;

  /* Interpolation_n-D: '<S254>/CLdeltafr' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_i[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_i[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_i[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_i[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_i[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_i[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_i[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_i[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_i[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_i[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_i[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_i[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_i[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_i[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_i[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_i[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_i[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_i[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_i[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_i[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_i[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_i[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_i[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_i[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_i[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_i[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_i[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_i[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_i[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_i[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_i[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_i[7] = 0;
  }

  mainV03_56_B.CLdeltafr = intrp8d_s32dla_pw(bpIndex_i, frac_i,
    mainV03_56_P.CLdeltafr_Table, mainV03_56_P.CLdeltafr_dimSize,
    mainV03_56_P.CLdeltafr_maxIndex);

  /* End of Interpolation_n-D: '<S254>/CLdeltafr' */

  /* Product: '<S254>/Product10' */
  mainV03_56_B.Product10 = mainV03_56_B.CLdeltafr *
    mainV03_56_B.ActuatorsCmd.deltafr;

  /* Interpolation_n-D: '<S254>/CLdeltafrl' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_j[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_j[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_j[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_j[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_j[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_j[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_j[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_j[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_j[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_j[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_j[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_j[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_j[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_j[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_j[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_j[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_j[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_j[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_j[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_j[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_j[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_j[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_j[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_j[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_j[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_j[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_j[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_j[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_j[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_j[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_j[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_j[7] = 0;
  }

  mainV03_56_B.CLdeltafrl = intrp8d_s32dla_pw(bpIndex_j, frac_j,
    mainV03_56_P.CLdeltafrl_Table, mainV03_56_P.CLdeltafrl_dimSize,
    mainV03_56_P.CLdeltafrl_maxIndex);

  /* End of Interpolation_n-D: '<S254>/CLdeltafrl' */

  /* Product: '<S254>/Product11' */
  mainV03_56_B.Product11 = mainV03_56_B.CLdeltafrl *
    mainV03_56_B.ActuatorsCmd.deltafl;

  /* Sum: '<S254>/Add1' */
  mainV03_56_B.Add1 = (((((mainV03_56_B.CL0 + mainV03_56_B.Product6_e) +
    mainV03_56_B.Product7_j) + mainV03_56_B.Product8_h) +
                        mainV03_56_B.Product9_j) + mainV03_56_B.Product10) +
    mainV03_56_B.Product11;

  /* SignalConversion: '<S242>/TmpSignal ConversionAtProductInport2' */
  mainV03_56_B.TmpSignalConversionAtProductInport2[0] = mainV03_56_B.Add;
  mainV03_56_B.TmpSignalConversionAtProductInport2[1] = mainV03_56_B.Add_g;
  mainV03_56_B.TmpSignalConversionAtProductInport2[2] = mainV03_56_B.Add1;

  /* Product: '<S242>/Product' incorporates:
   *  Math: '<S242>/Transpose'
   */
  for (rtb_Sum1_eh = 0; rtb_Sum1_eh < 3; rtb_Sum1_eh++) {
    mainV03_56_B.Product_i[rtb_Sum1_eh] = 0.0;
    mainV03_56_B.Product_i[rtb_Sum1_eh] += mainV03_56_B.VectorConcatenate_l[3 *
      rtb_Sum1_eh] * mainV03_56_B.TmpSignalConversionAtProductInport2[0];
    mainV03_56_B.Product_i[rtb_Sum1_eh] += mainV03_56_B.VectorConcatenate_l[3 *
      rtb_Sum1_eh + 1] * mainV03_56_B.TmpSignalConversionAtProductInport2[1];
    mainV03_56_B.Product_i[rtb_Sum1_eh] += mainV03_56_B.VectorConcatenate_l[3 *
      rtb_Sum1_eh + 2] * mainV03_56_B.TmpSignalConversionAtProductInport2[2];
  }

  /* End of Product: '<S242>/Product' */
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
  }

  /* Product: '<S287>/Product' */
  mainV03_56_B.Product_f = mainV03_56_B.Sum_ba[0] * mainV03_56_B.Sum_ba[0];

  /* Product: '<S287>/Product1' */
  mainV03_56_B.Product1_n = mainV03_56_B.Sum_ba[1] * mainV03_56_B.Sum_ba[1];

  /* Product: '<S287>/Product2' */
  mainV03_56_B.Product2_mn = mainV03_56_B.Sum_ba[2] * mainV03_56_B.Sum_ba[2];

  /* Sum: '<S287>/Sum' */
  mainV03_56_B.Sum_c4 = (mainV03_56_B.Product_f + mainV03_56_B.Product1_n) +
    mainV03_56_B.Product2_mn;

  /* Product: '<S260>/Product2' */
  mainV03_56_B.Product2_e = mainV03_56_B.Sum_c4 *
    mainV03_56_B.envData.AirDensity;

  /* Gain: '<S260>/1//2rhoV^2' */
  mainV03_56_B.u2rhoV2 = mainV03_56_P.u2rhoV2_Gain * mainV03_56_B.Product2_e;

  /* Gain: '<S252>/reference area' */
  mainV03_56_B.referencearea = mainV03_56_P.AerodynamicForcesandMoments_S *
    mainV03_56_B.u2rhoV2;

  /* Gain: '<S252>/coefAdjust' */
  mainV03_56_B.coefAdjust[0] = mainV03_56_P.coefAdjust_Gain[0] *
    mainV03_56_B.TmpSignalConversionAtProductInport2[0];

  /* Product: '<S252>/Product' */
  mainV03_56_B.Product_kd[0] = mainV03_56_B.referencearea *
    mainV03_56_B.coefAdjust[0];

  /* Gain: '<S252>/coefAdjust' */
  mainV03_56_B.coefAdjust[1] = mainV03_56_P.coefAdjust_Gain[1] *
    mainV03_56_B.TmpSignalConversionAtProductInport2[1];

  /* Product: '<S252>/Product' */
  mainV03_56_B.Product_kd[1] = mainV03_56_B.referencearea *
    mainV03_56_B.coefAdjust[1];

  /* Gain: '<S252>/coefAdjust' */
  mainV03_56_B.coefAdjust[2] = mainV03_56_P.coefAdjust_Gain[2] *
    mainV03_56_B.TmpSignalConversionAtProductInport2[2];

  /* Product: '<S252>/Product' */
  mainV03_56_B.Product_kd[2] = mainV03_56_B.referencearea *
    mainV03_56_B.coefAdjust[2];

  /* Trigonometry: '<S269>/Incidence' */
  rtb_Rn = rt_atan2d_snf(mainV03_56_B.Sum_ba[2], mainV03_56_B.Sum_ba[0]);

  /* Trigonometry: '<S268>/Trigonometric Function2' */
  rtb_UnaryMinus_j = sin(rtb_Rn);
  rtb_UnaryMinus_o = cos(rtb_Rn);

  /* SignalConversion: '<S270>/ConcatBufferAtVector ConcatenateIn1' */
  mainV03_56_B.VectorConcatenate_d[0] = rtb_UnaryMinus_o;
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Constant: '<S268>/Constant1' */
    mainV03_56_B.VectorConcatenate_d[1] = mainV03_56_P.Constant1_Value_f3;
  }

  /* UnaryMinus: '<S268>/Unary Minus' */
  mainV03_56_B.VectorConcatenate_d[2] = -rtb_UnaryMinus_j;
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Constant: '<S268>/Constant' */
    mainV03_56_B.VectorConcatenate_d[3] = mainV03_56_P.Constant_Value_bv;

    /* Constant: '<S268>/Constant2' */
    mainV03_56_B.VectorConcatenate_d[4] = mainV03_56_P.Constant2_Value_lrz;

    /* Constant: '<S268>/Constant4' */
    mainV03_56_B.VectorConcatenate_d[5] = mainV03_56_P.Constant4_Value_i;
  }

  /* SignalConversion: '<S270>/ConcatBufferAtVector ConcatenateIn7' */
  mainV03_56_B.VectorConcatenate_d[6] = rtb_UnaryMinus_j;
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Constant: '<S268>/Constant3' */
    mainV03_56_B.VectorConcatenate_d[7] = mainV03_56_P.Constant3_Value_o;
  }

  /* SignalConversion: '<S270>/ConcatBufferAtVector ConcatenateIn9' */
  mainV03_56_B.VectorConcatenate_d[8] = rtb_UnaryMinus_o;

  /* Sum: '<S252>/Sum' */
  mainV03_56_B.Sum_hf[0] = mainV03_56_B.plantData.CG[0] -
    mainV03_56_B.plantData.CG[0];
  mainV03_56_B.Sum_hf[1] = mainV03_56_B.plantData.CG[1] -
    mainV03_56_B.plantData.CG[1];
  mainV03_56_B.Sum_hf[2] = mainV03_56_B.plantData.CG[2] -
    mainV03_56_B.plantData.CG[2];

  /* Product: '<S263>/Product' */
  for (rtb_Sum1_eh = 0; rtb_Sum1_eh < 3; rtb_Sum1_eh++) {
    mainV03_56_B.Product_e[rtb_Sum1_eh] = 0.0;
    mainV03_56_B.Product_e[rtb_Sum1_eh] +=
      mainV03_56_B.VectorConcatenate_d[rtb_Sum1_eh] * mainV03_56_B.Sum_hf[0];
    mainV03_56_B.Product_e[rtb_Sum1_eh] +=
      mainV03_56_B.VectorConcatenate_d[rtb_Sum1_eh + 3] * mainV03_56_B.Sum_hf[1];
    mainV03_56_B.Product_e[rtb_Sum1_eh] +=
      mainV03_56_B.VectorConcatenate_d[rtb_Sum1_eh + 6] * mainV03_56_B.Sum_hf[2];
  }

  /* End of Product: '<S263>/Product' */

  /* Product: '<S266>/i x j' */
  mainV03_56_B.ixj_p = mainV03_56_B.Product_kd[0] * mainV03_56_B.Product_e[1];

  /* Product: '<S266>/j x k' */
  mainV03_56_B.jxk_do = mainV03_56_B.Product_kd[1] * mainV03_56_B.Product_e[2];

  /* Product: '<S266>/k x i' */
  mainV03_56_B.kxi_gg = mainV03_56_B.Product_kd[2] * mainV03_56_B.Product_e[0];

  /* Product: '<S267>/i x k' */
  mainV03_56_B.ixk_og = mainV03_56_B.Product_kd[0] * mainV03_56_B.Product_e[2];

  /* Product: '<S267>/j x i' */
  mainV03_56_B.jxi_c = mainV03_56_B.Product_kd[1] * mainV03_56_B.Product_e[0];

  /* Product: '<S267>/k x j' */
  mainV03_56_B.kxj_d = mainV03_56_B.Product_kd[2] * mainV03_56_B.Product_e[1];

  /* Sum: '<S262>/Sum' */
  mainV03_56_B.Sum_j[0] = mainV03_56_B.jxk_do - mainV03_56_B.kxj_d;
  mainV03_56_B.Sum_j[1] = mainV03_56_B.kxi_gg - mainV03_56_B.ixk_og;
  mainV03_56_B.Sum_j[2] = mainV03_56_B.ixj_p - mainV03_56_B.jxi_c;
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Constant: '<S274>/Constant' */
    mainV03_56_B.VectorConcatenate_d0[3] = mainV03_56_P.Constant_Value_e5;

    /* Constant: '<S274>/Constant1' */
    mainV03_56_B.VectorConcatenate_d0[1] = mainV03_56_P.Constant1_Value_lq;

    /* Constant: '<S274>/Constant2' */
    mainV03_56_B.VectorConcatenate_d0[4] = mainV03_56_P.Constant2_Value_nf;

    /* Constant: '<S274>/Constant3' */
    mainV03_56_B.VectorConcatenate_d0[7] = mainV03_56_P.Constant3_Value_mn;

    /* Constant: '<S274>/Constant4' */
    mainV03_56_B.VectorConcatenate_d0[5] = mainV03_56_P.Constant4_Value_pb;
  }

  /* Trigonometry: '<S275>/Incidence' */
  rtb_Rn = rt_atan2d_snf(mainV03_56_B.Sum_ba[2], mainV03_56_B.Sum_ba[0]);

  /* Trigonometry: '<S274>/Trigonometric Function2' */
  rtb_UnaryMinus_j = sin(rtb_Rn);
  rtb_UnaryMinus_o = cos(rtb_Rn);

  /* SignalConversion: '<S276>/ConcatBufferAtVector ConcatenateIn1' */
  mainV03_56_B.VectorConcatenate_d0[0] = rtb_UnaryMinus_o;

  /* SignalConversion: '<S276>/ConcatBufferAtVector ConcatenateIn7' */
  mainV03_56_B.VectorConcatenate_d0[6] = rtb_UnaryMinus_j;

  /* SignalConversion: '<S276>/ConcatBufferAtVector ConcatenateIn9' */
  mainV03_56_B.VectorConcatenate_d0[8] = rtb_UnaryMinus_o;

  /* UnaryMinus: '<S274>/Unary Minus' */
  mainV03_56_B.VectorConcatenate_d0[2] = -rtb_UnaryMinus_j;

  /* Math: '<S264>/Transpose' */
  for (rtb_Sum1_eh = 0; rtb_Sum1_eh < 3; rtb_Sum1_eh++) {
    mainV03_56_B.Transpose_h[3 * rtb_Sum1_eh] =
      mainV03_56_B.VectorConcatenate_d0[rtb_Sum1_eh];
    mainV03_56_B.Transpose_h[1 + 3 * rtb_Sum1_eh] =
      mainV03_56_B.VectorConcatenate_d0[rtb_Sum1_eh + 3];
    mainV03_56_B.Transpose_h[2 + 3 * rtb_Sum1_eh] =
      mainV03_56_B.VectorConcatenate_d0[rtb_Sum1_eh + 6];
  }

  /* End of Math: '<S264>/Transpose' */

  /* Product: '<S264>/Product' */
  for (rtb_Sum1_eh = 0; rtb_Sum1_eh < 3; rtb_Sum1_eh++) {
    mainV03_56_B.Product_g0[rtb_Sum1_eh] = 0.0;
    mainV03_56_B.Product_g0[rtb_Sum1_eh] += mainV03_56_B.Transpose_h[rtb_Sum1_eh]
      * mainV03_56_B.Product_kd[0];
    mainV03_56_B.Product_g0[rtb_Sum1_eh] += mainV03_56_B.Transpose_h[rtb_Sum1_eh
      + 3] * mainV03_56_B.Product_kd[1];
    mainV03_56_B.Product_g0[rtb_Sum1_eh] += mainV03_56_B.Transpose_h[rtb_Sum1_eh
      + 6] * mainV03_56_B.Product_kd[2];
  }

  /* End of Product: '<S264>/Product' */
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Constant: '<S280>/Constant' */
    mainV03_56_B.VectorConcatenate_j[3] = mainV03_56_P.Constant_Value_dh;

    /* Constant: '<S280>/Constant1' */
    mainV03_56_B.VectorConcatenate_j[1] = mainV03_56_P.Constant1_Value_j;

    /* Constant: '<S280>/Constant2' */
    mainV03_56_B.VectorConcatenate_j[4] = mainV03_56_P.Constant2_Value_ho;

    /* Constant: '<S280>/Constant3' */
    mainV03_56_B.VectorConcatenate_j[7] = mainV03_56_P.Constant3_Value_i;

    /* Constant: '<S280>/Constant4' */
    mainV03_56_B.VectorConcatenate_j[5] = mainV03_56_P.Constant4_Value_k;
  }

  /* Trigonometry: '<S281>/Incidence' */
  rtb_Rn = rt_atan2d_snf(mainV03_56_B.Sum_ba[2], mainV03_56_B.Sum_ba[0]);

  /* Trigonometry: '<S280>/Trigonometric Function2' */
  rtb_UnaryMinus_j = sin(rtb_Rn);
  rtb_UnaryMinus_o = cos(rtb_Rn);

  /* SignalConversion: '<S282>/ConcatBufferAtVector ConcatenateIn1' */
  mainV03_56_B.VectorConcatenate_j[0] = rtb_UnaryMinus_o;

  /* SignalConversion: '<S282>/ConcatBufferAtVector ConcatenateIn7' */
  mainV03_56_B.VectorConcatenate_j[6] = rtb_UnaryMinus_j;

  /* SignalConversion: '<S282>/ConcatBufferAtVector ConcatenateIn9' */
  mainV03_56_B.VectorConcatenate_j[8] = rtb_UnaryMinus_o;

  /* UnaryMinus: '<S280>/Unary Minus' */
  mainV03_56_B.VectorConcatenate_j[2] = -rtb_UnaryMinus_j;

  /* Math: '<S265>/Transpose' */
  for (rtb_Sum1_eh = 0; rtb_Sum1_eh < 3; rtb_Sum1_eh++) {
    mainV03_56_B.Transpose_e[3 * rtb_Sum1_eh] =
      mainV03_56_B.VectorConcatenate_j[rtb_Sum1_eh];
    mainV03_56_B.Transpose_e[1 + 3 * rtb_Sum1_eh] =
      mainV03_56_B.VectorConcatenate_j[rtb_Sum1_eh + 3];
    mainV03_56_B.Transpose_e[2 + 3 * rtb_Sum1_eh] =
      mainV03_56_B.VectorConcatenate_j[rtb_Sum1_eh + 6];
  }

  /* End of Math: '<S265>/Transpose' */

  /* Interpolation_n-D: '<S256>/Clbeta' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_k[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_k[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_k[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_k[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_k[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_k[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_k[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_k[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_k[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_k[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_k[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_k[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_k[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_k[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_k[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_k[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_k[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_k[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_k[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_k[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_k[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_k[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_k[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_k[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_k[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_k[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_k[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_k[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_k[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_k[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_k[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_k[7] = 0;
  }

  mainV03_56_B.Clbeta = intrp8d_s32dla_pw(bpIndex_k, frac_k,
    mainV03_56_P.Clbeta_Table, mainV03_56_P.Clbeta_dimSize,
    mainV03_56_P.Clbeta_maxIndex);

  /* End of Interpolation_n-D: '<S256>/Clbeta' */

  /* Product: '<S256>/Product6' */
  mainV03_56_B.Product6_h = mainV03_56_B.Clbeta * mainV03_56_B.plantData.beta;

  /* Interpolation_n-D: '<S256>/Clp' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_l[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_l[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_l[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_l[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_l[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_l[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_l[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_l[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_l[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_l[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_l[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_l[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_l[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_l[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_l[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_l[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_l[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_l[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_l[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_l[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_l[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_l[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_l[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_l[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_l[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_l[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_l[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_l[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_l[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_l[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_l[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_l[7] = 0;
  }

  mainV03_56_B.Clp = intrp8d_s32dla_pw(bpIndex_l, frac_l, mainV03_56_P.Clp_Table,
    mainV03_56_P.Clp_dimSize, mainV03_56_P.Clp_maxIndex);

  /* End of Interpolation_n-D: '<S256>/Clp' */

  /* Product: '<S256>/Product7' */
  mainV03_56_B.Product7_e = mainV03_56_B.Clp *
    mainV03_56_B.plantData.Omega_body.p;

  /* Interpolation_n-D: '<S256>/Clr' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_m[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_m[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_m[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_m[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_m[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_m[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_m[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_m[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_m[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_m[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_m[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_m[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_m[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_m[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_m[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_m[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_m[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_m[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_m[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_m[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_m[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_m[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_m[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_m[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_m[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_m[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_m[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_m[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_m[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_m[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_m[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_m[7] = 0;
  }

  mainV03_56_B.Clr = intrp8d_s32dla_pw(bpIndex_m, frac_m, mainV03_56_P.Clr_Table,
    mainV03_56_P.Clr_dimSize, mainV03_56_P.Clr_maxIndex);

  /* End of Interpolation_n-D: '<S256>/Clr' */

  /* Product: '<S256>/Product8' */
  mainV03_56_B.Product8_j = mainV03_56_B.Clr *
    mainV03_56_B.plantData.Omega_body.r;

  /* Interpolation_n-D: '<S256>/Cldeltar ' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_n[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_n[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_n[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_n[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_n[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_n[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_n[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_n[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_n[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_n[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_n[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_n[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_n[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_n[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_n[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_n[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_n[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_n[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_n[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_n[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_n[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_n[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_n[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_n[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_n[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_n[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_n[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_n[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_n[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_n[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_n[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_n[7] = 0;
  }

  mainV03_56_B.Cldeltar = intrp8d_s32dla_pw(bpIndex_n, frac_n,
    mainV03_56_P.Cldeltar_Table, mainV03_56_P.Cldeltar_dimSize,
    mainV03_56_P.Cldeltar_maxIndex);

  /* End of Interpolation_n-D: '<S256>/Cldeltar ' */

  /* Product: '<S256>/Product9' */
  mainV03_56_B.Product9_d = mainV03_56_B.Cldeltar *
    mainV03_56_B.ActuatorsCmd.deltar;

  /* Interpolation_n-D: '<S256>/Cldeltafr' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_o[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_o[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_o[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_o[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_o[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_o[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_o[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_o[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_o[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_o[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_o[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_o[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_o[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_o[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_o[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_o[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_o[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_o[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_o[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_o[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_o[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_o[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_o[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_o[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_o[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_o[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_o[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_o[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_o[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_o[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_o[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_o[7] = 0;
  }

  mainV03_56_B.Cldeltafr = intrp8d_s32dla_pw(bpIndex_o, frac_o,
    mainV03_56_P.Cldeltafr_Table, mainV03_56_P.Cldeltafr_dimSize,
    mainV03_56_P.Cldeltafr_maxIndex);

  /* End of Interpolation_n-D: '<S256>/Cldeltafr' */

  /* Product: '<S256>/Product10' */
  mainV03_56_B.Product10_d = mainV03_56_B.Cldeltafr *
    mainV03_56_B.ActuatorsCmd.deltafr;

  /* Interpolation_n-D: '<S256>/Cldeltafl' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_p[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_p[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_p[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_p[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_p[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_p[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_p[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_p[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_p[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_p[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_p[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_p[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_p[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_p[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_p[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_p[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_p[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_p[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_p[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_p[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_p[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_p[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_p[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_p[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_p[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_p[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_p[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_p[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_p[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_p[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_p[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_p[7] = 0;
  }

  mainV03_56_B.Cldeltafl = intrp8d_s32dla_pw(bpIndex_p, frac_p,
    mainV03_56_P.Cldeltafl_Table, mainV03_56_P.Cldeltafl_dimSize,
    mainV03_56_P.Cldeltafl_maxIndex);

  /* End of Interpolation_n-D: '<S256>/Cldeltafl' */

  /* Product: '<S256>/Product11' */
  mainV03_56_B.Product11_d = mainV03_56_B.Cldeltafl *
    mainV03_56_B.ActuatorsCmd.deltafl;

  /* Sum: '<S256>/Add1' */
  mainV03_56_B.Add1_f = ((((mainV03_56_B.Product6_h + mainV03_56_B.Product7_e) +
    mainV03_56_B.Product8_j) + mainV03_56_B.Product9_d) +
    mainV03_56_B.Product10_d) + mainV03_56_B.Product11_d;

  /* Interpolation_n-D: '<S257>/Cm0' */
  frac_q[0] = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  frac_q[1] = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  frac_q[2] = mainV03_56_B.analysisCasesBus_j.alt.falt;
  frac_q[3] = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  frac_q[4] = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  frac_q[5] = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  frac_q[6] = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  frac_q[7] = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 7) {
    bpIndex_q[0] = 7;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_q[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_q[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 0) {
    bpIndex_q[1] = 0;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_q[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_q[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 0) {
    bpIndex_q[2] = 0;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_q[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_q[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 0) {
    bpIndex_q[3] = 0;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_q[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_q[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 0) {
    bpIndex_q[4] = 0;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_q[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_q[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 0) {
    bpIndex_q[5] = 0;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_q[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_q[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 0) {
    bpIndex_q[6] = 0;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_q[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_q[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 0) {
    bpIndex_q[7] = 0;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_q[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_q[7] = 0;
  }

  mainV03_56_B.Cm0 = intrp8d_s32dl_pw(bpIndex_q, frac_q, mainV03_56_P.Cm0_Table,
    mainV03_56_P.Cm0_dimSize);

  /* End of Interpolation_n-D: '<S257>/Cm0' */

  /* Interpolation_n-D: '<S257>/Cmalpha' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_r[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_r[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_r[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_r[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_r[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_r[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_r[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_r[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_r[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_r[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_r[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_r[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_r[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_r[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_r[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_r[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_r[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_r[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_r[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_r[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_r[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_r[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_r[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_r[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_r[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_r[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_r[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_r[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_r[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_r[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_r[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_r[7] = 0;
  }

  mainV03_56_B.Cmalpha = intrp8d_s32dla_pw(bpIndex_r, frac_r,
    mainV03_56_P.Cmalpha_Table, mainV03_56_P.Cmalpha_dimSize,
    mainV03_56_P.Cmalpha_maxIndex);

  /* End of Interpolation_n-D: '<S257>/Cmalpha' */

  /* Product: '<S257>/Product6' */
  mainV03_56_B.Product6_ep = mainV03_56_B.Cmalpha * mainV03_56_B.plantData.alpha;

  /* Interpolation_n-D: '<S257>/Cmalpha_dot' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_s[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_s[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_s[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_s[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_s[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_s[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_s[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_s[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_s[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_s[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_s[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_s[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_s[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_s[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_s[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_s[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_s[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_s[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_s[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_s[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_s[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_s[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_s[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_s[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_s[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_s[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_s[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_s[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_s[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_s[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_s[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_s[7] = 0;
  }

  mainV03_56_B.Cmalpha_dot = intrp8d_s32dla_pw(bpIndex_s, frac_s,
    mainV03_56_P.Cmalpha_dot_Table, mainV03_56_P.Cmalpha_dot_dimSize,
    mainV03_56_P.Cmalpha_dot_maxIndex);

  /* End of Interpolation_n-D: '<S257>/Cmalpha_dot' */

  /* Product: '<S257>/Product7' */
  mainV03_56_B.Product7_k = mainV03_56_B.Cmalpha_dot *
    mainV03_56_B.plantData.alpha_dot;

  /* Interpolation_n-D: '<S257>/Cmq' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_t[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_t[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_t[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_t[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_t[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_t[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_t[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_t[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_t[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_t[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_t[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_t[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_t[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_t[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_t[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_t[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_t[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_t[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_t[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_t[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_t[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_t[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_t[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_t[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_t[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_t[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_t[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_t[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_t[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_t[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_t[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_t[7] = 0;
  }

  mainV03_56_B.Cmq = intrp8d_s32dla_pw(bpIndex_t, frac_t, mainV03_56_P.Cmq_Table,
    mainV03_56_P.Cmq_dimSize, mainV03_56_P.Cmq_maxIndex);

  /* End of Interpolation_n-D: '<S257>/Cmq' */

  /* Product: '<S257>/Product8' */
  mainV03_56_B.Product8_p = mainV03_56_B.Cmq *
    mainV03_56_B.plantData.Omega_body.q;

  /* Interpolation_n-D: '<S257>/Cmdeltae' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_u[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_u[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_u[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_u[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_u[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_u[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_u[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_u[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_u[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_u[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_u[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_u[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_u[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_u[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_u[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_u[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_u[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_u[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_u[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_u[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_u[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_u[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_u[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_u[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_u[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_u[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_u[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_u[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_u[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_u[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_u[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_u[7] = 0;
  }

  mainV03_56_B.Cmdeltae = intrp8d_s32dla_pw(bpIndex_u, frac_u,
    mainV03_56_P.Cmdeltae_Table, mainV03_56_P.Cmdeltae_dimSize,
    mainV03_56_P.Cmdeltae_maxIndex);

  /* End of Interpolation_n-D: '<S257>/Cmdeltae' */

  /* Product: '<S257>/Product9' */
  mainV03_56_B.Product9_p = mainV03_56_B.Cmdeltae *
    mainV03_56_B.ActuatorsCmd.deltae;

  /* Interpolation_n-D: '<S257>/Cmdeltafr' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_v[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_v[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_v[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_v[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_v[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_v[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_v[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_v[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_v[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_v[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_v[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_v[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_v[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_v[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_v[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_v[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_v[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_v[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_v[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_v[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_v[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_v[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_v[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_v[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_v[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_v[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_v[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_v[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_v[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_v[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_v[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_v[7] = 0;
  }

  mainV03_56_B.Cmdeltafr = intrp8d_s32dla_pw(bpIndex_v, frac_v,
    mainV03_56_P.Cmdeltafr_Table, mainV03_56_P.Cmdeltafr_dimSize,
    mainV03_56_P.Cmdeltafr_maxIndex);

  /* End of Interpolation_n-D: '<S257>/Cmdeltafr' */

  /* Product: '<S257>/Product10' */
  mainV03_56_B.Product10_j = mainV03_56_B.Cmdeltafr *
    mainV03_56_B.ActuatorsCmd.deltafr;

  /* Interpolation_n-D: '<S257>/Cmdeltafl' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_w[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_w[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_w[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_w[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_w[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_w[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_w[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_w[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_w[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_w[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_w[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_w[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_w[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_w[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_w[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_w[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_w[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_w[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_w[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_w[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_w[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_w[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_w[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_w[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_w[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_w[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_w[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_w[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_w[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_w[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_w[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_w[7] = 0;
  }

  mainV03_56_B.Cmdeltafl = intrp8d_s32dla_pw(bpIndex_w, frac_w,
    mainV03_56_P.Cmdeltafl_Table, mainV03_56_P.Cmdeltafl_dimSize,
    mainV03_56_P.Cmdeltafl_maxIndex);

  /* End of Interpolation_n-D: '<S257>/Cmdeltafl' */

  /* Product: '<S257>/Product11' */
  mainV03_56_B.Product11_a = mainV03_56_B.Cmdeltafl *
    mainV03_56_B.ActuatorsCmd.deltafl;

  /* Sum: '<S257>/Add2' */
  mainV03_56_B.Add2 = (((((mainV03_56_B.Cm0 + mainV03_56_B.Product6_ep) +
    mainV03_56_B.Product7_k) + mainV03_56_B.Product8_p) +
                        mainV03_56_B.Product9_p) + mainV03_56_B.Product10_j) +
    mainV03_56_B.Product11_a;

  /* Interpolation_n-D: '<S258>/Cnbeta' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_x[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_x[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_x[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_x[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_x[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_x[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_x[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_x[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_x[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_x[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_x[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_x[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_x[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_x[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_x[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_x[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_x[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_x[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_x[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_x[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_x[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_x[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_x[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_x[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_x[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_x[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_x[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_x[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_x[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_x[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_x[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_x[7] = 0;
  }

  mainV03_56_B.Cnbeta = intrp8d_s32dla_pw(bpIndex_x, frac_x,
    mainV03_56_P.Cnbeta_Table, mainV03_56_P.Cnbeta_dimSize,
    mainV03_56_P.Cnbeta_maxIndex);

  /* End of Interpolation_n-D: '<S258>/Cnbeta' */

  /* Product: '<S258>/Product6' */
  mainV03_56_B.Product6_m = mainV03_56_B.Cnbeta * mainV03_56_B.plantData.beta;

  /* Interpolation_n-D: '<S258>/Cnp' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_y[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_y[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_y[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_y[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_y[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_y[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_y[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_y[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_y[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_y[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_y[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_y[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_y[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_y[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_y[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_y[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_y[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_y[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_y[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_y[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_y[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_y[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_y[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_y[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_y[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_y[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_y[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_y[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_y[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_y[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_y[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_y[7] = 0;
  }

  mainV03_56_B.Cnp = intrp8d_s32dla_pw(bpIndex_y, frac_y, mainV03_56_P.Cnp_Table,
    mainV03_56_P.Cnp_dimSize, mainV03_56_P.Cnp_maxIndex);

  /* End of Interpolation_n-D: '<S258>/Cnp' */

  /* Product: '<S258>/Product7' */
  mainV03_56_B.Product7_f = mainV03_56_B.Cnp *
    mainV03_56_B.plantData.Omega_body.p;

  /* Interpolation_n-D: '<S258>/Cnr' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_z[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_z[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_z[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_z[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_z[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_z[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_z[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_z[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_z[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_z[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_z[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_z[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_z[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_z[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_z[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_z[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_z[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_z[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_z[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_z[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_z[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_z[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_z[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_z[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_z[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_z[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_z[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_z[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_z[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_z[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_z[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_z[7] = 0;
  }

  mainV03_56_B.Cnr = intrp8d_s32dla_pw(bpIndex_z, frac_z, mainV03_56_P.Cnr_Table,
    mainV03_56_P.Cnr_dimSize, mainV03_56_P.Cnr_maxIndex);

  /* End of Interpolation_n-D: '<S258>/Cnr' */

  /* Product: '<S258>/Product8' */
  mainV03_56_B.Product8_f = mainV03_56_B.Cnr *
    mainV03_56_B.plantData.Omega_body.r;

  /* Interpolation_n-D: '<S258>/Cndeltar ' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_10[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_10[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_10[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_10[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_10[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_10[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_10[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_10[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_10[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_10[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_10[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_10[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_10[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_10[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_10[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_10[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_10[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_10[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_10[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_10[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_10[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_10[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_10[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_10[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_10[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_10[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_10[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_10[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_10[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_10[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_10[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_10[7] = 0;
  }

  mainV03_56_B.Cndeltar = intrp8d_s32dla_pw(bpIndex_10, frac_10,
    mainV03_56_P.Cndeltar_Table, mainV03_56_P.Cndeltar_dimSize,
    mainV03_56_P.Cndeltar_maxIndex);

  /* End of Interpolation_n-D: '<S258>/Cndeltar ' */

  /* Product: '<S258>/Product9' */
  mainV03_56_B.Product9_e = mainV03_56_B.Cndeltar *
    mainV03_56_B.ActuatorsCmd.deltar;

  /* Interpolation_n-D: '<S258>/Cndeltafr' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_11[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_11[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_11[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_11[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_11[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_11[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_11[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_11[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_11[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_11[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_11[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_11[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_11[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_11[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_11[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_11[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_11[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_11[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_11[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_11[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_11[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_11[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_11[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_11[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_11[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_11[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_11[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_11[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_11[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_11[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_11[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_11[7] = 0;
  }

  mainV03_56_B.Cndeltafr = intrp8d_s32dla_pw(bpIndex_11, frac_11,
    mainV03_56_P.Cndeltafr_Table, mainV03_56_P.Cndeltafr_dimSize,
    mainV03_56_P.Cndeltafr_maxIndex);

  /* End of Interpolation_n-D: '<S258>/Cndeltafr' */

  /* Product: '<S258>/Product10' */
  mainV03_56_B.Product10_k = mainV03_56_B.Cndeltafr *
    mainV03_56_B.ActuatorsCmd.deltafr;

  /* Interpolation_n-D: '<S258>/Cndeltafl' */
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alpha.falpha;
  if (mainV03_56_B.analysisCasesBus_j.alpha.falpha < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alpha.falpha > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_12[0] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.beta.fbeta;
  if (mainV03_56_B.analysisCasesBus_j.beta.fbeta < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.beta.fbeta > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_12[1] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.alt.falt;
  if (mainV03_56_B.analysisCasesBus_j.alt.falt < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.alt.falt > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_12[2] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.xcg.fxcg;
  if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.xcg.fxcg > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_12[3] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltae.fdeltae;
  if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltae.fdeltae > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_12[4] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltar.fdeltar;
  if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltar.fdeltar > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_12[5] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr;
  if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafr.fdeltafr > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_12[6] = rtb_Rn;
  rtb_Rn = mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl;
  if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl < 0.0) {
    rtb_Rn = 0.0;
  } else {
    if (mainV03_56_B.analysisCasesBus_j.deltafl.fdeltafl > 1.0) {
      rtb_Rn = 1.0;
    }
  }

  frac_12[7] = rtb_Rn;
  if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha > 8) {
    bpIndex_12[0] = 8;
  } else if (mainV03_56_B.analysisCasesBus_j.alpha.idxalpha >= 0) {
    bpIndex_12[0] = mainV03_56_B.analysisCasesBus_j.alpha.idxalpha;
  } else {
    bpIndex_12[0] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta > 1) {
    bpIndex_12[1] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.beta.idxbeta >= 0) {
    bpIndex_12[1] = mainV03_56_B.analysisCasesBus_j.beta.idxbeta;
  } else {
    bpIndex_12[1] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.alt.idxalt > 1) {
    bpIndex_12[2] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.alt.idxalt >= 0) {
    bpIndex_12[2] = mainV03_56_B.analysisCasesBus_j.alt.idxalt;
  } else {
    bpIndex_12[2] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg > 1) {
    bpIndex_12[3] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.xcg.idxxcg >= 0) {
    bpIndex_12[3] = mainV03_56_B.analysisCasesBus_j.xcg.idxxcg;
  } else {
    bpIndex_12[3] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae > 1) {
    bpIndex_12[4] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae >= 0) {
    bpIndex_12[4] = mainV03_56_B.analysisCasesBus_j.deltae.idxdeltae;
  } else {
    bpIndex_12[4] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar > 1) {
    bpIndex_12[5] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar >= 0) {
    bpIndex_12[5] = mainV03_56_B.analysisCasesBus_j.deltar.idxdeltar;
  } else {
    bpIndex_12[5] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr > 1) {
    bpIndex_12[6] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr >= 0) {
    bpIndex_12[6] = mainV03_56_B.analysisCasesBus_j.deltafr.idxdeltafr;
  } else {
    bpIndex_12[6] = 0;
  }

  if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl > 1) {
    bpIndex_12[7] = 1;
  } else if (mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl >= 0) {
    bpIndex_12[7] = mainV03_56_B.analysisCasesBus_j.deltafl.idxdeltafl;
  } else {
    bpIndex_12[7] = 0;
  }

  mainV03_56_B.Cndeltafl = intrp8d_s32dla_pw(bpIndex_12, frac_12,
    mainV03_56_P.Cndeltafl_Table, mainV03_56_P.Cndeltafl_dimSize,
    mainV03_56_P.Cndeltafl_maxIndex);

  /* End of Interpolation_n-D: '<S258>/Cndeltafl' */

  /* Product: '<S258>/Product11' */
  mainV03_56_B.Product11_g = mainV03_56_B.Cndeltafl *
    mainV03_56_B.ActuatorsCmd.deltafl;

  /* Sum: '<S258>/Add1' */
  mainV03_56_B.Add1_e = ((((mainV03_56_B.Product6_m + mainV03_56_B.Product7_f) +
    mainV03_56_B.Product8_f) + mainV03_56_B.Product9_e) +
    mainV03_56_B.Product10_k) + mainV03_56_B.Product11_g;

  /* Product: '<S252>/Product1' incorporates:
   *  Constant: '<S252>/Constant'
   *  Constant: '<S252>/Constant1'
   */
  mainV03_56_B.Product1_k[0] = mainV03_56_P.AerodynamicForcesandMoments_b *
    mainV03_56_B.referencearea;
  mainV03_56_B.Product1_k[1] = mainV03_56_P.AerodynamicForcesandMoments_cbar *
    mainV03_56_B.referencearea;
  mainV03_56_B.Product1_k[2] = mainV03_56_P.AerodynamicForcesandMoments_b *
    mainV03_56_B.referencearea;

  /* Product: '<S252>/Product3' */
  mainV03_56_B.Product3_ix[0] = mainV03_56_B.Add1_f * mainV03_56_B.Product1_k[0];
  mainV03_56_B.Product3_ix[1] = mainV03_56_B.Add2 * mainV03_56_B.Product1_k[1];
  mainV03_56_B.Product3_ix[2] = mainV03_56_B.Add1_e * mainV03_56_B.Product1_k[2];

  /* Sum: '<S252>/Sum1' */
  mainV03_56_B.Sum1_g[0] = mainV03_56_B.Sum_j[0] + mainV03_56_B.Product3_ix[0];
  mainV03_56_B.Sum1_g[1] = mainV03_56_B.Sum_j[1] + mainV03_56_B.Product3_ix[1];
  mainV03_56_B.Sum1_g[2] = mainV03_56_B.Sum_j[2] + mainV03_56_B.Product3_ix[2];

  /* Product: '<S265>/Product' */
  for (rtb_Sum1_eh = 0; rtb_Sum1_eh < 3; rtb_Sum1_eh++) {
    mainV03_56_B.Product_fr[rtb_Sum1_eh] = 0.0;
    mainV03_56_B.Product_fr[rtb_Sum1_eh] += mainV03_56_B.Transpose_e[rtb_Sum1_eh]
      * mainV03_56_B.Sum1_g[0];
    mainV03_56_B.Product_fr[rtb_Sum1_eh] += mainV03_56_B.Transpose_e[rtb_Sum1_eh
      + 3] * mainV03_56_B.Sum1_g[1];
    mainV03_56_B.Product_fr[rtb_Sum1_eh] += mainV03_56_B.Transpose_e[rtb_Sum1_eh
      + 6] * mainV03_56_B.Sum1_g[2];
  }

  /* End of Product: '<S265>/Product' */

  /* InitialCondition: '<S242>/IC1' */
  if ((mainV03_56_DW.IC1_FirstOutputTime_k == (rtMinusInf)) ||
      (mainV03_56_DW.IC1_FirstOutputTime_k == mainV03_56_M->Timing.t[0])) {
    mainV03_56_DW.IC1_FirstOutputTime_k = mainV03_56_M->Timing.t[0];
    rtb_Shiftright[0] = mainV03_56_P.IC1_Value_d[0];
    rtb_Shiftright[1] = mainV03_56_P.IC1_Value_d[1];
    rtb_Shiftright[2] = mainV03_56_P.IC1_Value_d[2];
  } else {
    rtb_Shiftright[0] = mainV03_56_B.Product_g0[0];
    rtb_Shiftright[1] = mainV03_56_B.Product_g0[1];
    rtb_Shiftright[2] = mainV03_56_B.Product_g0[2];
  }

  /* End of InitialCondition: '<S242>/IC1' */

  /* InitialCondition: '<S242>/IC2' */
  if ((mainV03_56_DW.IC2_FirstOutputTime_g == (rtMinusInf)) ||
      (mainV03_56_DW.IC2_FirstOutputTime_g == mainV03_56_M->Timing.t[0])) {
    mainV03_56_DW.IC2_FirstOutputTime_g = mainV03_56_M->Timing.t[0];
    rtb_sincos_o2_p[0] = mainV03_56_P.IC2_Value_a[0];
    rtb_sincos_o2_p[1] = mainV03_56_P.IC2_Value_a[1];
    rtb_sincos_o2_p[2] = mainV03_56_P.IC2_Value_a[2];
  } else {
    rtb_sincos_o2_p[0] = mainV03_56_B.Product_fr[0];
    rtb_sincos_o2_p[1] = mainV03_56_B.Product_fr[1];
    rtb_sincos_o2_p[2] = mainV03_56_B.Product_fr[2];
  }

  /* End of InitialCondition: '<S242>/IC2' */

  /* Trigonometry: '<S300>/sincos' */
  rtb_sincos_o1_n_idx_0 = sin(mainV03_56_B.phithetapsi[0]);
  rtb_sincos_o2_i_idx_0 = cos(mainV03_56_B.phithetapsi[0]);
  rtb_sincos_o2_i_idx_1 = cos(mainV03_56_B.phithetapsi[1]);

  /* Fcn: '<S300>/phidot' incorporates:
   *  Trigonometry: '<S300>/sincos'
   */
  mainV03_56_B.phidot = (mainV03_56_B.pqr[1] * rtb_sincos_o1_n_idx_0 +
    mainV03_56_B.pqr[2] * rtb_sincos_o2_i_idx_0) * (sin
    (mainV03_56_B.phithetapsi[1]) / rtb_sincos_o2_i_idx_1) + mainV03_56_B.pqr[0];

  /* Fcn: '<S300>/psidot' */
  mainV03_56_B.psidot = (mainV03_56_B.pqr[1] * rtb_sincos_o1_n_idx_0 +
    mainV03_56_B.pqr[2] * rtb_sincos_o2_i_idx_0) / rtb_sincos_o2_i_idx_1;

  /* Fcn: '<S300>/thetadot' */
  mainV03_56_B.thetadot = mainV03_56_B.pqr[1] * rtb_sincos_o2_i_idx_0 -
    mainV03_56_B.pqr[2] * rtb_sincos_o1_n_idx_0;
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Constant: '<S246>/Ixx' */
    rtb_VectorConcatenate_i[0] = mainV03_56_P.Ixx_Value;

    /* SignalConversion: '<S341>/ConcatBufferAtVector ConcatenateIn2' incorporates:
     *  Constant: '<S246>/Ixy'
     *  UnaryMinus: '<S337>/Unary Minus2'
     */
    rtb_VectorConcatenate_i[1] = -mainV03_56_P.Ixy_Value;

    /* SignalConversion: '<S341>/ConcatBufferAtVector ConcatenateIn3' incorporates:
     *  Constant: '<S246>/Ixz'
     *  UnaryMinus: '<S337>/Unary Minus1'
     */
    rtb_VectorConcatenate_i[2] = -mainV03_56_P.Ixz_Value;

    /* SignalConversion: '<S341>/ConcatBufferAtVector ConcatenateIn4' incorporates:
     *  Constant: '<S246>/Ixy'
     *  UnaryMinus: '<S337>/Unary Minus2'
     */
    rtb_VectorConcatenate_i[3] = -mainV03_56_P.Ixy_Value;

    /* Constant: '<S246>/Iyy' */
    rtb_VectorConcatenate_i[4] = mainV03_56_P.Iyy_Value;

    /* SignalConversion: '<S341>/ConcatBufferAtVector ConcatenateIn6' incorporates:
     *  Constant: '<S246>/Iyz'
     *  UnaryMinus: '<S337>/Unary Minus'
     */
    rtb_VectorConcatenate_i[5] = -mainV03_56_P.Iyz_Value;

    /* SignalConversion: '<S341>/ConcatBufferAtVector ConcatenateIn7' incorporates:
     *  Constant: '<S246>/Ixz'
     *  UnaryMinus: '<S337>/Unary Minus1'
     */
    rtb_VectorConcatenate_i[6] = -mainV03_56_P.Ixz_Value;

    /* SignalConversion: '<S341>/ConcatBufferAtVector ConcatenateIn8' incorporates:
     *  Constant: '<S246>/Iyz'
     *  UnaryMinus: '<S337>/Unary Minus'
     */
    rtb_VectorConcatenate_i[7] = -mainV03_56_P.Iyz_Value;

    /* Constant: '<S246>/Izz' */
    rtb_VectorConcatenate_i[8] = mainV03_56_P.Izz_Value;

    /* Constant: '<S246>/dIxx//dt' */
    rtb_VectorConcatenate_n[0] = mainV03_56_P.dIxxdt_Value;

    /* SignalConversion: '<S342>/ConcatBufferAtVector ConcatenateIn2' incorporates:
     *  Constant: '<S246>/dIxy//dt'
     *  UnaryMinus: '<S338>/Unary Minus2'
     */
    rtb_VectorConcatenate_n[1] = -mainV03_56_P.dIxydt_Value;

    /* SignalConversion: '<S342>/ConcatBufferAtVector ConcatenateIn3' incorporates:
     *  Constant: '<S246>/dIxz//dt'
     *  UnaryMinus: '<S338>/Unary Minus1'
     */
    rtb_VectorConcatenate_n[2] = -mainV03_56_P.dIxzdt_Value;

    /* SignalConversion: '<S342>/ConcatBufferAtVector ConcatenateIn4' incorporates:
     *  Constant: '<S246>/dIxy//dt'
     *  UnaryMinus: '<S338>/Unary Minus2'
     */
    rtb_VectorConcatenate_n[3] = -mainV03_56_P.dIxydt_Value;

    /* Constant: '<S246>/dIyy//dt' */
    rtb_VectorConcatenate_n[4] = mainV03_56_P.dIyydt_Value;

    /* UnaryMinus: '<S338>/Unary Minus' incorporates:
     *  Constant: '<S246>/dIyz//dt'
     */
    rtb_Switch = -mainV03_56_P.dIyzdt_Value;

    /* SignalConversion: '<S342>/ConcatBufferAtVector ConcatenateIn6' incorporates:
     *  Constant: '<S246>/dIyz//dt'
     *  UnaryMinus: '<S338>/Unary Minus'
     */
    rtb_VectorConcatenate_n[5] = -mainV03_56_P.dIyzdt_Value;

    /* SignalConversion: '<S342>/ConcatBufferAtVector ConcatenateIn7' incorporates:
     *  Constant: '<S246>/dIxz//dt'
     *  UnaryMinus: '<S338>/Unary Minus1'
     */
    rtb_VectorConcatenate_n[6] = -mainV03_56_P.dIxzdt_Value;

    /* SignalConversion: '<S342>/ConcatBufferAtVector ConcatenateIn8' incorporates:
     *  Constant: '<S246>/dIyz//dt'
     *  UnaryMinus: '<S338>/Unary Minus'
     */
    rtb_VectorConcatenate_n[7] = -mainV03_56_P.dIyzdt_Value;

    /* Constant: '<S246>/dIzz//dt' */
    rtb_VectorConcatenate_n[8] = mainV03_56_P.dIzzdt_Value;
    for (rtb_Sum1_eh = 0; rtb_Sum1_eh < 3; rtb_Sum1_eh++) {
      /* Concatenate: '<S293>/Matrix Concatenation' */
      rtb_MatrixConcatenation[6 * rtb_Sum1_eh] = rtb_VectorConcatenate_i[3 *
        rtb_Sum1_eh];
      rtb_MatrixConcatenation[3 + 6 * rtb_Sum1_eh] = rtb_VectorConcatenate_n[3 *
        rtb_Sum1_eh];

      /* Selector: '<S292>/Selector' */
      mainV03_56_B.Selector_j[3 * rtb_Sum1_eh] = rtb_MatrixConcatenation[6 *
        rtb_Sum1_eh];

      /* Concatenate: '<S293>/Matrix Concatenation' */
      rtb_MatrixConcatenation[1 + 6 * rtb_Sum1_eh] = rtb_VectorConcatenate_i[3 *
        rtb_Sum1_eh + 1];
      rtb_MatrixConcatenation[4 + 6 * rtb_Sum1_eh] = rtb_VectorConcatenate_n[3 *
        rtb_Sum1_eh + 1];

      /* Selector: '<S292>/Selector' */
      mainV03_56_B.Selector_j[1 + 3 * rtb_Sum1_eh] = rtb_MatrixConcatenation[6 *
        rtb_Sum1_eh + 1];

      /* Concatenate: '<S293>/Matrix Concatenation' */
      rtb_MatrixConcatenation[2 + 6 * rtb_Sum1_eh] = rtb_VectorConcatenate_i[3 *
        rtb_Sum1_eh + 2];
      rtb_MatrixConcatenation[5 + 6 * rtb_Sum1_eh] = rtb_VectorConcatenate_n[3 *
        rtb_Sum1_eh + 2];

      /* Selector: '<S292>/Selector' */
      mainV03_56_B.Selector_j[2 + 3 * rtb_Sum1_eh] = rtb_MatrixConcatenation[6 *
        rtb_Sum1_eh + 2];
    }
  }

  /* Product: '<S302>/Product' */
  for (rtb_Sum1_eh = 0; rtb_Sum1_eh < 3; rtb_Sum1_eh++) {
    mainV03_56_B.Product_fi[rtb_Sum1_eh] = 0.0;
    mainV03_56_B.Product_fi[rtb_Sum1_eh] += mainV03_56_B.Selector_j[rtb_Sum1_eh]
      * mainV03_56_B.pqr[0];
    mainV03_56_B.Product_fi[rtb_Sum1_eh] += mainV03_56_B.Selector_j[rtb_Sum1_eh
      + 3] * mainV03_56_B.pqr[1];
    mainV03_56_B.Product_fi[rtb_Sum1_eh] += mainV03_56_B.Selector_j[rtb_Sum1_eh
      + 6] * mainV03_56_B.pqr[2];
  }

  /* End of Product: '<S302>/Product' */
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Selector: '<S292>/Selector1' */
    for (rtb_Sum1_eh = 0; rtb_Sum1_eh < 3; rtb_Sum1_eh++) {
      mainV03_56_B.Selector1_d[3 * rtb_Sum1_eh] = rtb_MatrixConcatenation[6 *
        rtb_Sum1_eh + 3];
      mainV03_56_B.Selector1_d[1 + 3 * rtb_Sum1_eh] = rtb_MatrixConcatenation[6 *
        rtb_Sum1_eh + 4];
      mainV03_56_B.Selector1_d[2 + 3 * rtb_Sum1_eh] = rtb_MatrixConcatenation[6 *
        rtb_Sum1_eh + 5];
    }

    /* End of Selector: '<S292>/Selector1' */
  }

  for (s155_iter = 0; s155_iter < 3; s155_iter++) {
    /* Product: '<S303>/Product' */
    mainV03_56_B.Product_ca[s155_iter] = 0.0;
    mainV03_56_B.Product_ca[s155_iter] += mainV03_56_B.Selector1_d[s155_iter] *
      mainV03_56_B.pqr[0];
    mainV03_56_B.Product_ca[s155_iter] += mainV03_56_B.Selector1_d[s155_iter + 3]
      * mainV03_56_B.pqr[1];
    mainV03_56_B.Product_ca[s155_iter] += mainV03_56_B.Selector1_d[s155_iter + 6]
      * mainV03_56_B.pqr[2];

    /* Sum: '<S349>/Sum1' incorporates:
     *  Constant: '<S349>/motor2 position '
     */
    mainV03_56_B.arm[s155_iter] = mainV03_56_P.motor2position_Value_k[s155_iter]
      - mainV03_56_B.plantData.CG[s155_iter];
  }

  /* Product: '<S355>/j x k' incorporates:
   *  Constant: '<S349>/Constant2'
   */
  mainV03_56_B.jxk_g = mainV03_56_B.arm[1] * mainV03_56_P.Constant2_Value_pm;
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Product: '<S247>/Product1' incorporates:
     *  Constant: '<S343>/Battery_Voltage'
     *  Constant: '<S343>/rendimiento_bateria'
     */
    rtb_Switch = mainV03_56_P.Battery_Voltage_Value *
      mainV03_56_P.rendimiento_bateria_Value;

    /* Product: '<S352>/Product' incorporates:
     *  Constant: '<S352>/eta_ESC'
     */
    mainV03_56_B.Product_p = rtb_Switch * mainV03_56_P.eta_ESC_Value;
  }

  /* Product: '<S352>/Product1' */
  mainV03_56_B.Product1_kx = mainV03_56_B.Product_p *
    mainV03_56_B.Throttle.Throttle1;

  /* PreLookup: '<S353>/Voltage Prelookup' */
  mainV03_56_B.VoltagePrelookup_o1 = plook_evenca(mainV03_56_B.Product1_kx,
    mainV03_56_P.VoltagePrelookup_BreakpointsData[0],
    mainV03_56_P.VoltagePrelookup_BreakpointsData[1] -
    mainV03_56_P.VoltagePrelookup_BreakpointsData[0], 149U,
    &mainV03_56_B.VoltagePrelookup_o2);

  /* PreLookup: '<S353>/Height Prelookup' */
  mainV03_56_B.HeightPrelookup_o1 = plook_evenca
    (mainV03_56_B.plantData.LLA.Altitude_m,
     mainV03_56_P.HeightPrelookup_BreakpointsData[0],
     mainV03_56_P.HeightPrelookup_BreakpointsData[1] -
     mainV03_56_P.HeightPrelookup_BreakpointsData[0], 4U,
     &mainV03_56_B.HeightPrelookup_o2);

  /* Sum: '<S349>/Sum' */
  mainV03_56_B.Sum_a5 = mainV03_56_B.plantData.V_body[0] -
    mainV03_56_B.envData.windVelocity[0];

  /* PreLookup: '<S353>/Vflight Prelookup' */
  mainV03_56_B.VflightPrelookup_o1 = plook_evenca(mainV03_56_B.Sum_a5,
    mainV03_56_P.VflightPrelookup_BreakpointsData_c[0],
    mainV03_56_P.VflightPrelookup_BreakpointsData_c[1] -
    mainV03_56_P.VflightPrelookup_BreakpointsData_c[0], 249U,
    &mainV03_56_B.VflightPrelookup_o2);

  /* Interpolation_n-D: '<S353>/Thrust Interpolation Using Prelookup' */
  frac_13[0] = mainV03_56_B.VoltagePrelookup_o2;
  frac_13[1] = mainV03_56_B.HeightPrelookup_o2;
  frac_13[2] = mainV03_56_B.VflightPrelookup_o2;
  bpIndex_13[0] = mainV03_56_B.VoltagePrelookup_o1;
  bpIndex_13[1] = mainV03_56_B.HeightPrelookup_o1;
  bpIndex_13[2] = mainV03_56_B.VflightPrelookup_o1;
  mainV03_56_B.Thrust = intrp3d_la_pw(bpIndex_13, frac_13,
    mainV03_56_P.ThrustInterpolationUsingPrelookup_Table,
    mainV03_56_P.ThrustInterpolationUsingPrelookup_dimSize,
    mainV03_56_P.ThrustInterpolationUsingPrelookup_maxIndex);

  /* Product: '<S355>/k x i' */
  mainV03_56_B.kxi_m = mainV03_56_B.arm[2] * mainV03_56_B.Thrust;

  /* Product: '<S355>/i x j' incorporates:
   *  Constant: '<S349>/Constant1'
   */
  mainV03_56_B.ixj_n = mainV03_56_B.arm[0] * mainV03_56_P.Constant1_Value_fx;

  /* Product: '<S356>/k x j' incorporates:
   *  Constant: '<S349>/Constant1'
   */
  mainV03_56_B.kxj_g = mainV03_56_B.arm[2] * mainV03_56_P.Constant1_Value_fx;

  /* Product: '<S356>/i x k' incorporates:
   *  Constant: '<S349>/Constant2'
   */
  mainV03_56_B.ixk_a = mainV03_56_B.arm[0] * mainV03_56_P.Constant2_Value_pm;

  /* Product: '<S356>/j x i' */
  mainV03_56_B.jxi_bv = mainV03_56_B.arm[1] * mainV03_56_B.Thrust;

  /* Sum: '<S351>/Sum' */
  mainV03_56_B.Sum_ac[0] = mainV03_56_B.jxk_g - mainV03_56_B.kxj_g;
  mainV03_56_B.Sum_ac[1] = mainV03_56_B.kxi_m - mainV03_56_B.ixk_a;
  mainV03_56_B.Sum_ac[2] = mainV03_56_B.ixj_n - mainV03_56_B.jxi_bv;

  /* Interpolation_n-D: '<S353>/RPM Interpolation Using Prelookup' */
  frac_14[0] = mainV03_56_B.VoltagePrelookup_o2;
  frac_14[1] = mainV03_56_B.HeightPrelookup_o2;
  frac_14[2] = mainV03_56_B.VflightPrelookup_o2;
  bpIndex_14[0] = mainV03_56_B.VoltagePrelookup_o1;
  bpIndex_14[1] = mainV03_56_B.HeightPrelookup_o1;
  bpIndex_14[2] = mainV03_56_B.VflightPrelookup_o1;
  mainV03_56_B.RPM = intrp3d_la_pw(bpIndex_14, frac_14,
    mainV03_56_P.RPMInterpolationUsingPrelookup_Table,
    mainV03_56_P.RPMInterpolationUsingPrelookup_dimSize,
    mainV03_56_P.RPMInterpolationUsingPrelookup_maxIndex);

  /* Sum: '<S354>/Sum' incorporates:
   *  Constant: '<S354>/Constant'
   */
  mainV03_56_B.Sum_kp = mainV03_56_B.RPM - mainV03_56_P.Constant_Value_nn;

  /* Switch: '<S354>/Switch' incorporates:
   *  Constant: '<S354>/Constant1'
   */
  if (mainV03_56_B.Sum_kp >= mainV03_56_P.Switch_Threshold_c) {
    /* PreLookup: '<S354>/Vflight Prelookup' */
    mainV03_56_B.VflightPrelookup_o1_f = plook_binca(mainV03_56_B.Sum_a5,
      mainV03_56_P.VflightPrelookup_BreakpointsData, 999U,
      &mainV03_56_B.VflightPrelookup_o2_dh);

    /* PreLookup: '<S354>/RPM Prelookup' */
    mainV03_56_B.RPMPrelookup_o1_f = plook_evenca(mainV03_56_B.RPM,
      mainV03_56_P.RPMPrelookup_BreakpointsData[0],
      mainV03_56_P.RPMPrelookup_BreakpointsData[1] -
      mainV03_56_P.RPMPrelookup_BreakpointsData[0], 14U,
      &mainV03_56_B.RPMPrelookup_o2_n);

    /* Interpolation_n-D: '<S354>/Torque Interpolation Using Prelookup' */
    frac_1i[0] = mainV03_56_B.RPMPrelookup_o2_n;
    frac_1i[1] = mainV03_56_B.VflightPrelookup_o2_dh;
    bpIndex_1i[0] = mainV03_56_B.RPMPrelookup_o1_f;
    bpIndex_1i[1] = mainV03_56_B.VflightPrelookup_o1_f;
    mainV03_56_B.Torque_az = intrp2d_la_pw(bpIndex_1i, frac_1i,
      mainV03_56_P.TorqueInterpolationUsingPrelookup_Table, 15U,
      mainV03_56_P.TorqueInterpolationUsingPrelookup_maxIndex);
    mainV03_56_B.Torque = mainV03_56_B.Torque_az;
  } else {
    mainV03_56_B.Torque = mainV03_56_P.Constant1_Value_jd;
  }

  /* End of Switch: '<S354>/Switch' */

  /* Sum: '<S349>/Sum2' incorporates:
   *  Constant: '<S349>/Constant3'
   *  Constant: '<S349>/Constant4'
   */
  mainV03_56_B.Sum2[0] = mainV03_56_B.Sum_ac[0] + mainV03_56_B.Torque;
  mainV03_56_B.Sum2[1] = mainV03_56_B.Sum_ac[1] +
    mainV03_56_P.Constant3_Value_nr;
  mainV03_56_B.Sum2[2] = mainV03_56_B.Sum_ac[2] +
    mainV03_56_P.Constant4_Value_ie;

  /* Sum: '<S363>/Sum1' incorporates:
   *  Constant: '<S363>/motor2 position '
   */
  mainV03_56_B.arm_h[0] = mainV03_56_P.motor2position_Value_i[0] -
    mainV03_56_B.plantData.CG[0];
  mainV03_56_B.arm_h[1] = mainV03_56_P.motor2position_Value_i[1] -
    mainV03_56_B.plantData.CG[1];
  mainV03_56_B.arm_h[2] = mainV03_56_P.motor2position_Value_i[2] -
    mainV03_56_B.plantData.CG[2];
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Product: '<S366>/Product' incorporates:
     *  Constant: '<S366>/eta_ESC'
     */
    mainV03_56_B.Product_dq = rtb_Switch * mainV03_56_P.eta_ESC_Value_d;
  }

  /* Product: '<S366>/Product1' */
  mainV03_56_B.Product1_c = mainV03_56_B.Product_dq *
    mainV03_56_B.Throttle.Throttle2;

  /* PreLookup: '<S367>/Voltage Prelookup' */
  mainV03_56_B.VoltagePrelookup_o1_a = plook_evenca(mainV03_56_B.Product1_c,
    mainV03_56_P.VoltagePrelookup_BreakpointsData_k[0],
    mainV03_56_P.VoltagePrelookup_BreakpointsData_k[1] -
    mainV03_56_P.VoltagePrelookup_BreakpointsData_k[0], 149U,
    &mainV03_56_B.VoltagePrelookup_o2_j);

  /* PreLookup: '<S367>/Height Prelookup' */
  mainV03_56_B.HeightPrelookup_o1_g = plook_evenca
    (mainV03_56_B.plantData.LLA.Altitude_m,
     mainV03_56_P.HeightPrelookup_BreakpointsData_c[0],
     mainV03_56_P.HeightPrelookup_BreakpointsData_c[1] -
     mainV03_56_P.HeightPrelookup_BreakpointsData_c[0], 4U,
     &mainV03_56_B.HeightPrelookup_o2_k);

  /* Sum: '<S363>/Sum' */
  mainV03_56_B.Sum_eo = mainV03_56_B.plantData.V_body[2] -
    mainV03_56_B.envData.windVelocity[2];

  /* PreLookup: '<S367>/Vflight Prelookup' */
  mainV03_56_B.VflightPrelookup_o1_e = plook_evenca(mainV03_56_B.Sum_eo,
    mainV03_56_P.VflightPrelookup_BreakpointsData_p[0],
    mainV03_56_P.VflightPrelookup_BreakpointsData_p[1] -
    mainV03_56_P.VflightPrelookup_BreakpointsData_p[0], 249U,
    &mainV03_56_B.VflightPrelookup_o2_d);

  /* Interpolation_n-D: '<S367>/Thrust Interpolation Using Prelookup' */
  frac_15[0] = mainV03_56_B.VoltagePrelookup_o2_j;
  frac_15[1] = mainV03_56_B.HeightPrelookup_o2_k;
  frac_15[2] = mainV03_56_B.VflightPrelookup_o2_d;
  bpIndex_15[0] = mainV03_56_B.VoltagePrelookup_o1_a;
  bpIndex_15[1] = mainV03_56_B.HeightPrelookup_o1_g;
  bpIndex_15[2] = mainV03_56_B.VflightPrelookup_o1_e;
  mainV03_56_B.Thrust_d = intrp3d_la_pw(bpIndex_15, frac_15,
    mainV03_56_P.ThrustInterpolationUsingPrelookup_Table_a,
    mainV03_56_P.ThrustInterpolationUsingPrelookup_dimSize_c,
    mainV03_56_P.ThrustInterpolationUsingPrelookup_maxIndex_g);

  /* UnaryMinus: '<S363>/Unary Minus' */
  rtb_UnaryMinus_j = -mainV03_56_B.Thrust_d;

  /* Product: '<S369>/j x k' */
  mainV03_56_B.jxk_o = mainV03_56_B.arm_h[1] * rtb_UnaryMinus_j;

  /* Product: '<S369>/k x i' incorporates:
   *  Constant: '<S363>/Constant2'
   */
  mainV03_56_B.kxi_gm = mainV03_56_B.arm_h[2] * mainV03_56_P.Constant2_Value_j;

  /* Product: '<S369>/i x j' incorporates:
   *  Constant: '<S363>/Constant1'
   */
  mainV03_56_B.ixj_c = mainV03_56_B.arm_h[0] * mainV03_56_P.Constant1_Value_c;

  /* Product: '<S370>/k x j' incorporates:
   *  Constant: '<S363>/Constant1'
   */
  mainV03_56_B.kxj_gt = mainV03_56_B.arm_h[2] * mainV03_56_P.Constant1_Value_c;

  /* Product: '<S370>/i x k' */
  mainV03_56_B.ixk_b = mainV03_56_B.arm_h[0] * rtb_UnaryMinus_j;

  /* Product: '<S370>/j x i' incorporates:
   *  Constant: '<S363>/Constant2'
   */
  mainV03_56_B.jxi_gj = mainV03_56_B.arm_h[1] * mainV03_56_P.Constant2_Value_j;

  /* Sum: '<S365>/Sum' */
  mainV03_56_B.Sum_h5[0] = mainV03_56_B.jxk_o - mainV03_56_B.kxj_gt;
  mainV03_56_B.Sum_h5[1] = mainV03_56_B.kxi_gm - mainV03_56_B.ixk_b;
  mainV03_56_B.Sum_h5[2] = mainV03_56_B.ixj_c - mainV03_56_B.jxi_gj;

  /* Interpolation_n-D: '<S367>/RPM Interpolation Using Prelookup' */
  frac_16[0] = mainV03_56_B.VoltagePrelookup_o2_j;
  frac_16[1] = mainV03_56_B.HeightPrelookup_o2_k;
  frac_16[2] = mainV03_56_B.VflightPrelookup_o2_d;
  bpIndex_16[0] = mainV03_56_B.VoltagePrelookup_o1_a;
  bpIndex_16[1] = mainV03_56_B.HeightPrelookup_o1_g;
  bpIndex_16[2] = mainV03_56_B.VflightPrelookup_o1_e;
  mainV03_56_B.RPM_g = intrp3d_la_pw(bpIndex_16, frac_16,
    mainV03_56_P.RPMInterpolationUsingPrelookup_Table_i,
    mainV03_56_P.RPMInterpolationUsingPrelookup_dimSize_g,
    mainV03_56_P.RPMInterpolationUsingPrelookup_maxIndex_m);

  /* Sum: '<S368>/Sum' incorporates:
   *  Constant: '<S368>/Constant'
   */
  mainV03_56_B.Sum_m4 = mainV03_56_B.RPM_g - mainV03_56_P.Constant_Value_b3;

  /* Switch: '<S368>/Switch' incorporates:
   *  Constant: '<S368>/Constant1'
   */
  if (mainV03_56_B.Sum_m4 >= mainV03_56_P.Switch_Threshold_b) {
    /* PreLookup: '<S368>/Vflight Prelookup' */
    mainV03_56_B.VflightPrelookup_o1_j = plook_binca(mainV03_56_B.Sum_eo,
      mainV03_56_P.VflightPrelookup_BreakpointsData_j, 999U,
      &mainV03_56_B.VflightPrelookup_o2_j);

    /* PreLookup: '<S368>/RPM Prelookup' */
    mainV03_56_B.RPMPrelookup_o1_cy = plook_evenca(mainV03_56_B.RPM_g,
      mainV03_56_P.RPMPrelookup_BreakpointsData_n[0],
      mainV03_56_P.RPMPrelookup_BreakpointsData_n[1] -
      mainV03_56_P.RPMPrelookup_BreakpointsData_n[0], 11U,
      &mainV03_56_B.RPMPrelookup_o2_ao);

    /* Interpolation_n-D: '<S368>/Torque Interpolation Using Prelookup' */
    frac_1j[0] = mainV03_56_B.RPMPrelookup_o2_ao;
    frac_1j[1] = mainV03_56_B.VflightPrelookup_o2_j;
    bpIndex_1j[0] = mainV03_56_B.RPMPrelookup_o1_cy;
    bpIndex_1j[1] = mainV03_56_B.VflightPrelookup_o1_j;
    mainV03_56_B.Torque_a = intrp2d_la_pw(bpIndex_1j, frac_1j,
      mainV03_56_P.TorqueInterpolationUsingPrelookup_Table_b, 12U,
      mainV03_56_P.TorqueInterpolationUsingPrelookup_maxIndex_e);
    mainV03_56_B.Torque_h = mainV03_56_B.Torque_a;
  } else {
    mainV03_56_B.Torque_h = mainV03_56_P.Constant1_Value_lt;
  }

  /* End of Switch: '<S368>/Switch' */

  /* Sum: '<S363>/Sum2' incorporates:
   *  Constant: '<S363>/Constant3'
   *  Constant: '<S363>/Constant4'
   */
  mainV03_56_B.Sum2_a[0] = mainV03_56_B.Sum_h5[0] +
    mainV03_56_P.Constant4_Value_j;
  mainV03_56_B.Sum2_a[1] = mainV03_56_B.Sum_h5[1] +
    mainV03_56_P.Constant3_Value_j;
  mainV03_56_B.Sum2_a[2] = mainV03_56_B.Sum_h5[2] + mainV03_56_B.Torque_h;

  /* Sum: '<S377>/Sum1' incorporates:
   *  Constant: '<S377>/motor2 position '
   */
  mainV03_56_B.arm_e[0] = mainV03_56_P.motor2position_Value_n[0] -
    mainV03_56_B.plantData.CG[0];
  mainV03_56_B.arm_e[1] = mainV03_56_P.motor2position_Value_n[1] -
    mainV03_56_B.plantData.CG[1];
  mainV03_56_B.arm_e[2] = mainV03_56_P.motor2position_Value_n[2] -
    mainV03_56_B.plantData.CG[2];
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Product: '<S380>/Product' incorporates:
     *  Constant: '<S380>/eta_ESC'
     */
    mainV03_56_B.Product_j = rtb_Switch * mainV03_56_P.eta_ESC_Value_n;
  }

  /* Product: '<S380>/Product1' */
  mainV03_56_B.Product1_e = mainV03_56_B.Product_j *
    mainV03_56_B.Throttle.Throttle3;

  /* PreLookup: '<S381>/Voltage Prelookup' */
  mainV03_56_B.VoltagePrelookup_o1_l = plook_evenca(mainV03_56_B.Product1_e,
    mainV03_56_P.VoltagePrelookup_BreakpointsData_f[0],
    mainV03_56_P.VoltagePrelookup_BreakpointsData_f[1] -
    mainV03_56_P.VoltagePrelookup_BreakpointsData_f[0], 149U,
    &mainV03_56_B.VoltagePrelookup_o2_a);

  /* PreLookup: '<S381>/Height Prelookup' */
  mainV03_56_B.HeightPrelookup_o1_i = plook_evenca
    (mainV03_56_B.plantData.LLA.Altitude_m,
     mainV03_56_P.HeightPrelookup_BreakpointsData_b[0],
     mainV03_56_P.HeightPrelookup_BreakpointsData_b[1] -
     mainV03_56_P.HeightPrelookup_BreakpointsData_b[0], 4U,
     &mainV03_56_B.HeightPrelookup_o2_j);

  /* Sum: '<S377>/Sum' */
  mainV03_56_B.Sum_nq = mainV03_56_B.plantData.V_body[2] -
    mainV03_56_B.envData.windVelocity[2];

  /* PreLookup: '<S381>/Vflight Prelookup' */
  mainV03_56_B.VflightPrelookup_o1_m = plook_evenca(mainV03_56_B.Sum_nq,
    mainV03_56_P.VflightPrelookup_BreakpointsData_e[0],
    mainV03_56_P.VflightPrelookup_BreakpointsData_e[1] -
    mainV03_56_P.VflightPrelookup_BreakpointsData_e[0], 249U,
    &mainV03_56_B.VflightPrelookup_o2_f);

  /* Interpolation_n-D: '<S381>/Thrust Interpolation Using Prelookup' */
  frac_17[0] = mainV03_56_B.VoltagePrelookup_o2_a;
  frac_17[1] = mainV03_56_B.HeightPrelookup_o2_j;
  frac_17[2] = mainV03_56_B.VflightPrelookup_o2_f;
  bpIndex_17[0] = mainV03_56_B.VoltagePrelookup_o1_l;
  bpIndex_17[1] = mainV03_56_B.HeightPrelookup_o1_i;
  bpIndex_17[2] = mainV03_56_B.VflightPrelookup_o1_m;
  mainV03_56_B.Thrust_l = intrp3d_la_pw(bpIndex_17, frac_17,
    mainV03_56_P.ThrustInterpolationUsingPrelookup_Table_k,
    mainV03_56_P.ThrustInterpolationUsingPrelookup_dimSize_b,
    mainV03_56_P.ThrustInterpolationUsingPrelookup_maxIndex_o);

  /* UnaryMinus: '<S377>/Unary Minus' */
  rtb_UnaryMinus_o = -mainV03_56_B.Thrust_l;

  /* Product: '<S383>/j x k' */
  mainV03_56_B.jxk_b2 = mainV03_56_B.arm_e[1] * rtb_UnaryMinus_o;

  /* Product: '<S383>/k x i' incorporates:
   *  Constant: '<S377>/Constant2'
   */
  mainV03_56_B.kxi_go = mainV03_56_B.arm_e[2] * mainV03_56_P.Constant2_Value_d4;

  /* Product: '<S383>/i x j' incorporates:
   *  Constant: '<S377>/Constant1'
   */
  mainV03_56_B.ixj_dn = mainV03_56_B.arm_e[0] * mainV03_56_P.Constant1_Value_cp;

  /* Product: '<S384>/k x j' incorporates:
   *  Constant: '<S377>/Constant1'
   */
  mainV03_56_B.kxj_e = mainV03_56_B.arm_e[2] * mainV03_56_P.Constant1_Value_cp;

  /* Product: '<S384>/i x k' */
  mainV03_56_B.ixk_op = mainV03_56_B.arm_e[0] * rtb_UnaryMinus_o;

  /* Product: '<S384>/j x i' incorporates:
   *  Constant: '<S377>/Constant2'
   */
  mainV03_56_B.jxi_d = mainV03_56_B.arm_e[1] * mainV03_56_P.Constant2_Value_d4;

  /* Sum: '<S379>/Sum' */
  mainV03_56_B.Sum_j1[0] = mainV03_56_B.jxk_b2 - mainV03_56_B.kxj_e;
  mainV03_56_B.Sum_j1[1] = mainV03_56_B.kxi_go - mainV03_56_B.ixk_op;
  mainV03_56_B.Sum_j1[2] = mainV03_56_B.ixj_dn - mainV03_56_B.jxi_d;

  /* Interpolation_n-D: '<S381>/RPM Interpolation Using Prelookup' */
  frac_18[0] = mainV03_56_B.VoltagePrelookup_o2_a;
  frac_18[1] = mainV03_56_B.HeightPrelookup_o2_j;
  frac_18[2] = mainV03_56_B.VflightPrelookup_o2_f;
  bpIndex_18[0] = mainV03_56_B.VoltagePrelookup_o1_l;
  bpIndex_18[1] = mainV03_56_B.HeightPrelookup_o1_i;
  bpIndex_18[2] = mainV03_56_B.VflightPrelookup_o1_m;
  mainV03_56_B.RPM_f = intrp3d_la_pw(bpIndex_18, frac_18,
    mainV03_56_P.RPMInterpolationUsingPrelookup_Table_f,
    mainV03_56_P.RPMInterpolationUsingPrelookup_dimSize_f,
    mainV03_56_P.RPMInterpolationUsingPrelookup_maxIndex_mm);

  /* Sum: '<S382>/Sum' incorporates:
   *  Constant: '<S382>/Constant'
   */
  mainV03_56_B.Sum_cl = mainV03_56_B.RPM_f - mainV03_56_P.Constant_Value_hf;

  /* Switch: '<S382>/Switch' incorporates:
   *  Constant: '<S382>/Constant1'
   */
  if (mainV03_56_B.Sum_cl >= mainV03_56_P.Switch_Threshold_g) {
    /* PreLookup: '<S382>/Vflight Prelookup' */
    mainV03_56_B.VflightPrelookup_o1_gx = plook_binca(mainV03_56_B.Sum_nq,
      mainV03_56_P.VflightPrelookup_BreakpointsData_i, 999U,
      &mainV03_56_B.VflightPrelookup_o2_a);

    /* PreLookup: '<S382>/RPM Prelookup' */
    mainV03_56_B.RPMPrelookup_o1_c = plook_evenca(mainV03_56_B.RPM_f,
      mainV03_56_P.RPMPrelookup_BreakpointsData_f[0],
      mainV03_56_P.RPMPrelookup_BreakpointsData_f[1] -
      mainV03_56_P.RPMPrelookup_BreakpointsData_f[0], 11U,
      &mainV03_56_B.RPMPrelookup_o2_a);

    /* Interpolation_n-D: '<S382>/Torque Interpolation Using Prelookup' */
    frac_1k[0] = mainV03_56_B.RPMPrelookup_o2_a;
    frac_1k[1] = mainV03_56_B.VflightPrelookup_o2_a;
    bpIndex_1k[0] = mainV03_56_B.RPMPrelookup_o1_c;
    bpIndex_1k[1] = mainV03_56_B.VflightPrelookup_o1_gx;
    mainV03_56_B.Torque_hw = intrp2d_la_pw(bpIndex_1k, frac_1k,
      mainV03_56_P.TorqueInterpolationUsingPrelookup_Table_i, 12U,
      mainV03_56_P.TorqueInterpolationUsingPrelookup_maxIndex_b);
    mainV03_56_B.Torque_n = mainV03_56_B.Torque_hw;
  } else {
    mainV03_56_B.Torque_n = mainV03_56_P.Constant1_Value_ln;
  }

  /* End of Switch: '<S382>/Switch' */

  /* Sum: '<S377>/Sum2' incorporates:
   *  Constant: '<S377>/Constant3'
   *  Constant: '<S377>/Constant4'
   *  UnaryMinus: '<S377>/Unary Minus1'
   */
  mainV03_56_B.Sum2_h[0] = mainV03_56_B.Sum_j1[0] +
    mainV03_56_P.Constant4_Value_o;
  mainV03_56_B.Sum2_h[1] = mainV03_56_B.Sum_j1[1] +
    mainV03_56_P.Constant3_Value_pv;
  mainV03_56_B.Sum2_h[2] = mainV03_56_B.Sum_j1[2] + -mainV03_56_B.Torque_n;

  /* Sum: '<S391>/Sum1' incorporates:
   *  Constant: '<S391>/motor2 position '
   */
  mainV03_56_B.arm_p[0] = mainV03_56_P.motor2position_Value_nw[0] -
    mainV03_56_B.plantData.CG[0];
  mainV03_56_B.arm_p[1] = mainV03_56_P.motor2position_Value_nw[1] -
    mainV03_56_B.plantData.CG[1];
  mainV03_56_B.arm_p[2] = mainV03_56_P.motor2position_Value_nw[2] -
    mainV03_56_B.plantData.CG[2];
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Product: '<S394>/Product' incorporates:
     *  Constant: '<S394>/eta_ESC'
     */
    mainV03_56_B.Product_ch = rtb_Switch * mainV03_56_P.eta_ESC_Value_m;
  }

  /* Product: '<S394>/Product1' */
  mainV03_56_B.Product1_kd = mainV03_56_B.Product_ch *
    mainV03_56_B.Throttle.Throttle4;

  /* PreLookup: '<S395>/Voltage Prelookup' */
  mainV03_56_B.VoltagePrelookup_o1_i = plook_evenca(mainV03_56_B.Product1_kd,
    mainV03_56_P.VoltagePrelookup_BreakpointsData_m[0],
    mainV03_56_P.VoltagePrelookup_BreakpointsData_m[1] -
    mainV03_56_P.VoltagePrelookup_BreakpointsData_m[0], 149U,
    &mainV03_56_B.VoltagePrelookup_o2_k);

  /* PreLookup: '<S395>/Height Prelookup' */
  mainV03_56_B.HeightPrelookup_o1_d = plook_evenca
    (mainV03_56_B.plantData.LLA.Altitude_m,
     mainV03_56_P.HeightPrelookup_BreakpointsData_o[0],
     mainV03_56_P.HeightPrelookup_BreakpointsData_o[1] -
     mainV03_56_P.HeightPrelookup_BreakpointsData_o[0], 4U,
     &mainV03_56_B.HeightPrelookup_o2_e);

  /* Sum: '<S391>/Sum' */
  mainV03_56_B.Sum_fv = mainV03_56_B.plantData.V_body[2] -
    mainV03_56_B.envData.windVelocity[2];

  /* PreLookup: '<S395>/Vflight Prelookup' */
  mainV03_56_B.VflightPrelookup_o1_a = plook_evenca(mainV03_56_B.Sum_fv,
    mainV03_56_P.VflightPrelookup_BreakpointsData_jk[0],
    mainV03_56_P.VflightPrelookup_BreakpointsData_jk[1] -
    mainV03_56_P.VflightPrelookup_BreakpointsData_jk[0], 249U,
    &mainV03_56_B.VflightPrelookup_o2_fy);

  /* Interpolation_n-D: '<S395>/Thrust Interpolation Using Prelookup' */
  frac_19[0] = mainV03_56_B.VoltagePrelookup_o2_k;
  frac_19[1] = mainV03_56_B.HeightPrelookup_o2_e;
  frac_19[2] = mainV03_56_B.VflightPrelookup_o2_fy;
  bpIndex_19[0] = mainV03_56_B.VoltagePrelookup_o1_i;
  bpIndex_19[1] = mainV03_56_B.HeightPrelookup_o1_d;
  bpIndex_19[2] = mainV03_56_B.VflightPrelookup_o1_a;
  mainV03_56_B.Thrust_e = intrp3d_la_pw(bpIndex_19, frac_19,
    mainV03_56_P.ThrustInterpolationUsingPrelookup_Table_i,
    mainV03_56_P.ThrustInterpolationUsingPrelookup_dimSize_bl,
    mainV03_56_P.ThrustInterpolationUsingPrelookup_maxIndex_a);

  /* Product: '<S397>/j x k' incorporates:
   *  UnaryMinus: '<S391>/Unary Minus'
   */
  mainV03_56_B.jxk_a = mainV03_56_B.arm_p[1] * -mainV03_56_B.Thrust_e;

  /* Product: '<S397>/k x i' incorporates:
   *  Constant: '<S391>/Constant2'
   */
  mainV03_56_B.kxi_m5 = mainV03_56_B.arm_p[2] * mainV03_56_P.Constant2_Value_l5;

  /* Product: '<S397>/i x j' incorporates:
   *  Constant: '<S391>/Constant1'
   */
  mainV03_56_B.ixj_n0 = mainV03_56_B.arm_p[0] * mainV03_56_P.Constant1_Value_ea;

  /* Product: '<S398>/k x j' incorporates:
   *  Constant: '<S391>/Constant1'
   */
  mainV03_56_B.kxj_p = mainV03_56_B.arm_p[2] * mainV03_56_P.Constant1_Value_ea;

  /* Product: '<S398>/i x k' incorporates:
   *  UnaryMinus: '<S391>/Unary Minus'
   */
  mainV03_56_B.ixk_h1 = mainV03_56_B.arm_p[0] * -mainV03_56_B.Thrust_e;

  /* Product: '<S398>/j x i' incorporates:
   *  Constant: '<S391>/Constant2'
   */
  mainV03_56_B.jxi_m = mainV03_56_B.arm_p[1] * mainV03_56_P.Constant2_Value_l5;

  /* Sum: '<S393>/Sum' */
  mainV03_56_B.Sum_l[0] = mainV03_56_B.jxk_a - mainV03_56_B.kxj_p;
  mainV03_56_B.Sum_l[1] = mainV03_56_B.kxi_m5 - mainV03_56_B.ixk_h1;
  mainV03_56_B.Sum_l[2] = mainV03_56_B.ixj_n0 - mainV03_56_B.jxi_m;

  /* Interpolation_n-D: '<S395>/RPM Interpolation Using Prelookup' */
  frac_1a[0] = mainV03_56_B.VoltagePrelookup_o2_k;
  frac_1a[1] = mainV03_56_B.HeightPrelookup_o2_e;
  frac_1a[2] = mainV03_56_B.VflightPrelookup_o2_fy;
  bpIndex_1a[0] = mainV03_56_B.VoltagePrelookup_o1_i;
  bpIndex_1a[1] = mainV03_56_B.HeightPrelookup_o1_d;
  bpIndex_1a[2] = mainV03_56_B.VflightPrelookup_o1_a;
  mainV03_56_B.RPM_i = intrp3d_la_pw(bpIndex_1a, frac_1a,
    mainV03_56_P.RPMInterpolationUsingPrelookup_Table_p,
    mainV03_56_P.RPMInterpolationUsingPrelookup_dimSize_o,
    mainV03_56_P.RPMInterpolationUsingPrelookup_maxIndex_c);

  /* Sum: '<S396>/Sum' incorporates:
   *  Constant: '<S396>/Constant'
   */
  mainV03_56_B.Sum_pt = mainV03_56_B.RPM_i - mainV03_56_P.Constant_Value_i2;

  /* Switch: '<S396>/Switch' incorporates:
   *  Constant: '<S396>/Constant1'
   */
  if (mainV03_56_B.Sum_pt >= mainV03_56_P.Switch_Threshold_a) {
    /* PreLookup: '<S396>/Vflight Prelookup' */
    mainV03_56_B.VflightPrelookup_o1_h = plook_binca(mainV03_56_B.Sum_fv,
      mainV03_56_P.VflightPrelookup_BreakpointsData_b, 999U,
      &mainV03_56_B.VflightPrelookup_o2_k);

    /* PreLookup: '<S396>/RPM Prelookup' */
    mainV03_56_B.RPMPrelookup_o1_k = plook_evenca(mainV03_56_B.RPM_i,
      mainV03_56_P.RPMPrelookup_BreakpointsData_m[0],
      mainV03_56_P.RPMPrelookup_BreakpointsData_m[1] -
      mainV03_56_P.RPMPrelookup_BreakpointsData_m[0], 11U,
      &mainV03_56_B.RPMPrelookup_o2_e);

    /* Interpolation_n-D: '<S396>/Torque Interpolation Using Prelookup' */
    frac_1l[0] = mainV03_56_B.RPMPrelookup_o2_e;
    frac_1l[1] = mainV03_56_B.VflightPrelookup_o2_k;
    bpIndex_1l[0] = mainV03_56_B.RPMPrelookup_o1_k;
    bpIndex_1l[1] = mainV03_56_B.VflightPrelookup_o1_h;
    mainV03_56_B.Torque_pl = intrp2d_la_pw(bpIndex_1l, frac_1l,
      mainV03_56_P.TorqueInterpolationUsingPrelookup_Table_h, 12U,
      mainV03_56_P.TorqueInterpolationUsingPrelookup_maxIndex_g);
    mainV03_56_B.Torque_c = mainV03_56_B.Torque_pl;
  } else {
    mainV03_56_B.Torque_c = mainV03_56_P.Constant1_Value_ct;
  }

  /* End of Switch: '<S396>/Switch' */

  /* Sum: '<S391>/Sum2' incorporates:
   *  Constant: '<S391>/Constant3'
   *  Constant: '<S391>/Constant4'
   */
  mainV03_56_B.Sum2_m[0] = mainV03_56_B.Sum_l[0] +
    mainV03_56_P.Constant4_Value_dl;
  mainV03_56_B.Sum2_m[1] = mainV03_56_B.Sum_l[1] +
    mainV03_56_P.Constant3_Value_nz;
  mainV03_56_B.Sum2_m[2] = mainV03_56_B.Sum_l[2] + mainV03_56_B.Torque_c;

  /* Sum: '<S405>/Sum1' incorporates:
   *  Constant: '<S405>/motor2 position '
   */
  mainV03_56_B.arm_a[0] = mainV03_56_P.motor2position_Value_l[0] -
    mainV03_56_B.plantData.CG[0];
  mainV03_56_B.arm_a[1] = mainV03_56_P.motor2position_Value_l[1] -
    mainV03_56_B.plantData.CG[1];
  mainV03_56_B.arm_a[2] = mainV03_56_P.motor2position_Value_l[2] -
    mainV03_56_B.plantData.CG[2];
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Product: '<S408>/Product' incorporates:
     *  Constant: '<S408>/eta_ESC'
     */
    mainV03_56_B.Product_h = rtb_Switch * mainV03_56_P.eta_ESC_Value_mv;
  }

  /* Product: '<S408>/Product1' */
  mainV03_56_B.Product1_ef = mainV03_56_B.Product_h *
    mainV03_56_B.Throttle.Throttle5;

  /* PreLookup: '<S409>/Voltage Prelookup' */
  mainV03_56_B.VoltagePrelookup_o1_a2 = plook_evenca(mainV03_56_B.Product1_ef,
    mainV03_56_P.VoltagePrelookup_BreakpointsData_d[0],
    mainV03_56_P.VoltagePrelookup_BreakpointsData_d[1] -
    mainV03_56_P.VoltagePrelookup_BreakpointsData_d[0], 149U,
    &mainV03_56_B.VoltagePrelookup_o2_f);

  /* PreLookup: '<S409>/Height Prelookup' */
  mainV03_56_B.HeightPrelookup_o1_l = plook_evenca
    (mainV03_56_B.plantData.LLA.Altitude_m,
     mainV03_56_P.HeightPrelookup_BreakpointsData_h[0],
     mainV03_56_P.HeightPrelookup_BreakpointsData_h[1] -
     mainV03_56_P.HeightPrelookup_BreakpointsData_h[0], 4U,
     &mainV03_56_B.HeightPrelookup_o2_o);

  /* Sum: '<S405>/Sum' */
  mainV03_56_B.Sum_ar = mainV03_56_B.plantData.V_body[2] -
    mainV03_56_B.envData.windVelocity[2];

  /* PreLookup: '<S409>/Vflight Prelookup' */
  mainV03_56_B.VflightPrelookup_o1_g = plook_evenca(mainV03_56_B.Sum_ar,
    mainV03_56_P.VflightPrelookup_BreakpointsData_jw[0],
    mainV03_56_P.VflightPrelookup_BreakpointsData_jw[1] -
    mainV03_56_P.VflightPrelookup_BreakpointsData_jw[0], 249U,
    &mainV03_56_B.VflightPrelookup_o2_m);

  /* Interpolation_n-D: '<S409>/Thrust Interpolation Using Prelookup' */
  frac_1b[0] = mainV03_56_B.VoltagePrelookup_o2_f;
  frac_1b[1] = mainV03_56_B.HeightPrelookup_o2_o;
  frac_1b[2] = mainV03_56_B.VflightPrelookup_o2_m;
  bpIndex_1b[0] = mainV03_56_B.VoltagePrelookup_o1_a2;
  bpIndex_1b[1] = mainV03_56_B.HeightPrelookup_o1_l;
  bpIndex_1b[2] = mainV03_56_B.VflightPrelookup_o1_g;
  mainV03_56_B.Thrust_m = intrp3d_la_pw(bpIndex_1b, frac_1b,
    mainV03_56_P.ThrustInterpolationUsingPrelookup_Table_m,
    mainV03_56_P.ThrustInterpolationUsingPrelookup_dimSize_n,
    mainV03_56_P.ThrustInterpolationUsingPrelookup_maxIndex_p);

  /* Product: '<S411>/j x k' incorporates:
   *  UnaryMinus: '<S405>/Unary Minus'
   */
  mainV03_56_B.jxk_e = mainV03_56_B.arm_a[1] * -mainV03_56_B.Thrust_m;

  /* Product: '<S411>/k x i' incorporates:
   *  Constant: '<S405>/Constant2'
   */
  mainV03_56_B.kxi_fl = mainV03_56_B.arm_a[2] * mainV03_56_P.Constant2_Value_lh;

  /* Product: '<S411>/i x j' incorporates:
   *  Constant: '<S405>/Constant1'
   */
  mainV03_56_B.ixj_eu = mainV03_56_B.arm_a[0] * mainV03_56_P.Constant1_Value_np;

  /* Product: '<S412>/k x j' incorporates:
   *  Constant: '<S405>/Constant1'
   */
  mainV03_56_B.kxj_m5 = mainV03_56_B.arm_a[2] * mainV03_56_P.Constant1_Value_np;

  /* Product: '<S412>/i x k' incorporates:
   *  UnaryMinus: '<S405>/Unary Minus'
   */
  mainV03_56_B.ixk_nc = mainV03_56_B.arm_a[0] * -mainV03_56_B.Thrust_m;

  /* Product: '<S412>/j x i' incorporates:
   *  Constant: '<S405>/Constant2'
   */
  mainV03_56_B.jxi_b4 = mainV03_56_B.arm_a[1] * mainV03_56_P.Constant2_Value_lh;

  /* Sum: '<S407>/Sum' */
  mainV03_56_B.Sum_ab[0] = mainV03_56_B.jxk_e - mainV03_56_B.kxj_m5;
  mainV03_56_B.Sum_ab[1] = mainV03_56_B.kxi_fl - mainV03_56_B.ixk_nc;
  mainV03_56_B.Sum_ab[2] = mainV03_56_B.ixj_eu - mainV03_56_B.jxi_b4;

  /* Interpolation_n-D: '<S409>/RPM Interpolation Using Prelookup' */
  frac_1c[0] = mainV03_56_B.VoltagePrelookup_o2_f;
  frac_1c[1] = mainV03_56_B.HeightPrelookup_o2_o;
  frac_1c[2] = mainV03_56_B.VflightPrelookup_o2_m;
  bpIndex_1c[0] = mainV03_56_B.VoltagePrelookup_o1_a2;
  bpIndex_1c[1] = mainV03_56_B.HeightPrelookup_o1_l;
  bpIndex_1c[2] = mainV03_56_B.VflightPrelookup_o1_g;
  mainV03_56_B.RPM_k = intrp3d_la_pw(bpIndex_1c, frac_1c,
    mainV03_56_P.RPMInterpolationUsingPrelookup_Table_k,
    mainV03_56_P.RPMInterpolationUsingPrelookup_dimSize_o5,
    mainV03_56_P.RPMInterpolationUsingPrelookup_maxIndex_h);

  /* Sum: '<S410>/Sum' incorporates:
   *  Constant: '<S410>/Constant'
   */
  mainV03_56_B.Sum_b0 = mainV03_56_B.RPM_k - mainV03_56_P.Constant_Value_h4;

  /* Switch: '<S410>/Switch' incorporates:
   *  Constant: '<S410>/Constant1'
   */
  if (mainV03_56_B.Sum_b0 >= mainV03_56_P.Switch_Threshold_g3) {
    /* PreLookup: '<S410>/Vflight Prelookup' */
    mainV03_56_B.VflightPrelookup_o1_d = plook_binca(mainV03_56_B.Sum_ar,
      mainV03_56_P.VflightPrelookup_BreakpointsData_m, 999U,
      &mainV03_56_B.VflightPrelookup_o2_n);

    /* PreLookup: '<S410>/RPM Prelookup' */
    mainV03_56_B.RPMPrelookup_o1 = plook_evenca(mainV03_56_B.RPM_k,
      mainV03_56_P.RPMPrelookup_BreakpointsData_ms[0],
      mainV03_56_P.RPMPrelookup_BreakpointsData_ms[1] -
      mainV03_56_P.RPMPrelookup_BreakpointsData_ms[0], 11U,
      &mainV03_56_B.RPMPrelookup_o2);

    /* Interpolation_n-D: '<S410>/Torque Interpolation Using Prelookup' */
    frac_1m[0] = mainV03_56_B.RPMPrelookup_o2;
    frac_1m[1] = mainV03_56_B.VflightPrelookup_o2_n;
    bpIndex_1m[0] = mainV03_56_B.RPMPrelookup_o1;
    bpIndex_1m[1] = mainV03_56_B.VflightPrelookup_o1_d;
    mainV03_56_B.Torque_b = intrp2d_la_pw(bpIndex_1m, frac_1m,
      mainV03_56_P.TorqueInterpolationUsingPrelookup_Table_p, 12U,
      mainV03_56_P.TorqueInterpolationUsingPrelookup_maxIndex_c);
    mainV03_56_B.Torque_p = mainV03_56_B.Torque_b;
  } else {
    mainV03_56_B.Torque_p = mainV03_56_P.Constant1_Value_nj;
  }

  /* End of Switch: '<S410>/Switch' */

  /* Sum: '<S405>/Sum2' incorporates:
   *  Constant: '<S405>/Constant3'
   *  Constant: '<S405>/Constant4'
   *  UnaryMinus: '<S405>/Unary Minus1'
   */
  mainV03_56_B.Sum2_c[0] = mainV03_56_B.Sum_ab[0] +
    mainV03_56_P.Constant4_Value_kq;
  mainV03_56_B.Sum2_c[1] = mainV03_56_B.Sum_ab[1] +
    mainV03_56_P.Constant3_Value_ki;
  mainV03_56_B.Sum2_c[2] = mainV03_56_B.Sum_ab[2] + -mainV03_56_B.Torque_p;

  /* Sum: '<S247>/Sum3' */
  mainV03_56_B.total_moment[0] = (((mainV03_56_B.Sum2[0] + mainV03_56_B.Sum2_a[0])
    + mainV03_56_B.Sum2_h[0]) + mainV03_56_B.Sum2_m[0]) + mainV03_56_B.Sum2_c[0];

  /* Sum: '<S4>/Sum1' */
  mainV03_56_B.Sum1_p3[0] = (rtb_sincos_o2_p[0] + mainV03_56_B.envData.Mground[0])
    + mainV03_56_B.total_moment[0];

  /* Sum: '<S247>/Sum3' */
  mainV03_56_B.total_moment[1] = (((mainV03_56_B.Sum2[1] + mainV03_56_B.Sum2_a[1])
    + mainV03_56_B.Sum2_h[1]) + mainV03_56_B.Sum2_m[1]) + mainV03_56_B.Sum2_c[1];

  /* Sum: '<S4>/Sum1' */
  mainV03_56_B.Sum1_p3[1] = (rtb_sincos_o2_p[1] + mainV03_56_B.envData.Mground[1])
    + mainV03_56_B.total_moment[1];

  /* Sum: '<S247>/Sum3' */
  mainV03_56_B.total_moment[2] = (((mainV03_56_B.Sum2[2] + mainV03_56_B.Sum2_a[2])
    + mainV03_56_B.Sum2_h[2]) + mainV03_56_B.Sum2_m[2]) + mainV03_56_B.Sum2_c[2];

  /* Sum: '<S4>/Sum1' */
  mainV03_56_B.Sum1_p3[2] = (rtb_sincos_o2_p[2] + mainV03_56_B.envData.Mground[2])
    + mainV03_56_B.total_moment[2];

  /* InitialCondition: '<S4>/IC9' */
  if ((mainV03_56_DW.IC9_FirstOutputTime == (rtMinusInf)) ||
      (mainV03_56_DW.IC9_FirstOutputTime == mainV03_56_M->Timing.t[0])) {
    mainV03_56_DW.IC9_FirstOutputTime = mainV03_56_M->Timing.t[0];
    rtb_sincos_o2_p[0] = mainV03_56_P.IC9_Value[0];
    rtb_sincos_o2_p[1] = mainV03_56_P.IC9_Value[1];
    rtb_sincos_o2_p[2] = mainV03_56_P.IC9_Value[2];
  } else {
    rtb_sincos_o2_p[0] = mainV03_56_B.Sum1_p3[0];
    rtb_sincos_o2_p[1] = mainV03_56_B.Sum1_p3[1];
    rtb_sincos_o2_p[2] = mainV03_56_B.Sum1_p3[2];
  }

  /* End of InitialCondition: '<S4>/IC9' */

  /* Product: '<S305>/j x k' */
  mainV03_56_B.jxk_cb = mainV03_56_B.pqr[1] * mainV03_56_B.Product_fi[2];

  /* Product: '<S305>/k x i' */
  mainV03_56_B.kxi_j = mainV03_56_B.pqr[2] * mainV03_56_B.Product_fi[0];

  /* Product: '<S305>/i x j' */
  mainV03_56_B.ixj_ey = mainV03_56_B.pqr[0] * mainV03_56_B.Product_fi[1];

  /* Product: '<S306>/k x j' */
  mainV03_56_B.kxj_d5 = mainV03_56_B.pqr[2] * mainV03_56_B.Product_fi[1];

  /* Product: '<S306>/i x k' */
  mainV03_56_B.ixk_fm = mainV03_56_B.pqr[0] * mainV03_56_B.Product_fi[2];

  /* Product: '<S306>/j x i' */
  mainV03_56_B.jxi_o = mainV03_56_B.pqr[1] * mainV03_56_B.Product_fi[0];

  /* Sum: '<S304>/Sum' */
  mainV03_56_B.Sum_bh[0] = mainV03_56_B.jxk_cb - mainV03_56_B.kxj_d5;
  mainV03_56_B.Sum_bh[1] = mainV03_56_B.kxi_j - mainV03_56_B.ixk_fm;
  mainV03_56_B.Sum_bh[2] = mainV03_56_B.ixj_ey - mainV03_56_B.jxi_o;

  /* Sum: '<S292>/Sum2' */
  mainV03_56_B.Sum2_mg[0] = (rtb_sincos_o2_p[0] - mainV03_56_B.Product_ca[0]) -
    mainV03_56_B.Sum_bh[0];
  mainV03_56_B.Sum2_mg[1] = (rtb_sincos_o2_p[1] - mainV03_56_B.Product_ca[1]) -
    mainV03_56_B.Sum_bh[1];
  mainV03_56_B.Sum2_mg[2] = (rtb_sincos_o2_p[2] - mainV03_56_B.Product_ca[2]) -
    mainV03_56_B.Sum_bh[2];
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Selector: '<S292>/Selector2' */
    for (rtb_Sum1_eh = 0; rtb_Sum1_eh < 3; rtb_Sum1_eh++) {
      mainV03_56_B.Selector2_n[3 * rtb_Sum1_eh] = rtb_MatrixConcatenation[6 *
        rtb_Sum1_eh];
      mainV03_56_B.Selector2_n[1 + 3 * rtb_Sum1_eh] = rtb_MatrixConcatenation[6 *
        rtb_Sum1_eh + 1];
      mainV03_56_B.Selector2_n[2 + 3 * rtb_Sum1_eh] = rtb_MatrixConcatenation[6 *
        rtb_Sum1_eh + 2];
    }

    /* End of Selector: '<S292>/Selector2' */
  }

  /* Product: '<S292>/Product2' */
  rt_mrdivide_U1d1x3_U2d3x3_Yd1x3_snf(mainV03_56_B.Sum2_mg,
    mainV03_56_B.Selector2_n, mainV03_56_B.Product2_ex);
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* SignalConversion: '<S307>/TmpSignal ConversionAtFor Each SubsystemInport2' incorporates:
     *  Constant: '<S246>/Vre_x (m//s)'
     *  Constant: '<S246>/Vre_y (m//s)'
     *  Constant: '<S246>/Vre_z (m//s)'
     */
    rtb_WhiteNoise_o[0] = mainV03_56_P.Vre_xms_Value;
    rtb_WhiteNoise_o[1] = mainV03_56_P.Vre_yms_Value;
    rtb_WhiteNoise_o[2] = mainV03_56_P.Vre_zms_Value;

    /* Outputs for Iterator SubSystem: '<S307>/For Each Subsystem' incorporates:
     *  ForEach: '<S308>/For Each'
     */
    for (ForEach_itr = 0; ForEach_itr < 1; ForEach_itr++) {
      /* ForEachSliceSelector: '<S308>/ImpSel_InsertedFor_Vre_at_outport_0' */
      rtb_Sum1_eh = ForEach_itr * 3;

      /* ForEachSliceAssignment: '<S308>/ImpAsg_InsertedFor_F_at_inport_0' incorporates:
       *  Constant: '<S246>/dmass//dt (kg//s)'
       *  ForEachSliceSelector: '<S308>/ImpSel_InsertedFor_Vre_at_outport_0'
       *  Product: '<S308>/Product'
       */
      rtb_ImpAsg_InsertedFor_F_at_inport_0_idx_0 = mainV03_56_P.dmassdtkgs_Value
        * rtb_WhiteNoise_o[rtb_Sum1_eh];
      rtb_ImpAsg_InsertedFor_F_at_inport_0_idx_1 = rtb_WhiteNoise_o[1 +
        rtb_Sum1_eh] * mainV03_56_P.dmassdtkgs_Value;
      rtb_ImpAsg_InsertedFor_F_at_inport_0_idx_2 = rtb_WhiteNoise_o[2 +
        rtb_Sum1_eh] * mainV03_56_P.dmassdtkgs_Value;
    }

    /* End of Outputs for SubSystem: '<S307>/For Each Subsystem' */

    /* Sum: '<S307>/Sum of Elements' */
    mainV03_56_B.SumofElements[0] = rtb_ImpAsg_InsertedFor_F_at_inport_0_idx_0;
    mainV03_56_B.SumofElements[1] = rtb_ImpAsg_InsertedFor_F_at_inport_0_idx_1;
    mainV03_56_B.SumofElements[2] = rtb_ImpAsg_InsertedFor_F_at_inport_0_idx_2;
  }

  /* InitialCondition: '<S4>/IC13' */
  if ((mainV03_56_DW.IC13_FirstOutputTime == (rtMinusInf)) ||
      (mainV03_56_DW.IC13_FirstOutputTime == mainV03_56_M->Timing.t[0])) {
    mainV03_56_DW.IC13_FirstOutputTime = mainV03_56_M->Timing.t[0];
    rtb_Shiftright[0] = mainV03_56_P.IC13_Value[0];
    rtb_Shiftright[1] = mainV03_56_P.IC13_Value[1];
    rtb_Shiftright[2] = mainV03_56_P.IC13_Value[2];
  }

  /* End of InitialCondition: '<S4>/IC13' */

  /* Sum: '<S247>/Sum2' incorporates:
   *  Constant: '<S349>/Constant1'
   *  Constant: '<S349>/Constant2'
   *  Constant: '<S363>/Constant1'
   *  Constant: '<S363>/Constant2'
   *  Constant: '<S377>/Constant1'
   *  Constant: '<S377>/Constant2'
   *  Constant: '<S391>/Constant1'
   *  Constant: '<S391>/Constant2'
   *  Constant: '<S405>/Constant1'
   *  Constant: '<S405>/Constant2'
   *  UnaryMinus: '<S391>/Unary Minus'
   *  UnaryMinus: '<S405>/Unary Minus'
   */
  mainV03_56_B.total_thrust[0] = (((mainV03_56_B.Thrust +
    mainV03_56_P.Constant2_Value_j) + mainV03_56_P.Constant2_Value_d4) +
    mainV03_56_P.Constant2_Value_l5) + mainV03_56_P.Constant2_Value_lh;
  mainV03_56_B.total_thrust[1] = (((mainV03_56_P.Constant1_Value_fx +
    mainV03_56_P.Constant1_Value_c) + mainV03_56_P.Constant1_Value_cp) +
    mainV03_56_P.Constant1_Value_ea) + mainV03_56_P.Constant1_Value_np;
  mainV03_56_B.total_thrust[2] = (((mainV03_56_P.Constant2_Value_pm +
    rtb_UnaryMinus_j) + rtb_UnaryMinus_o) + -mainV03_56_B.Thrust_e) +
    -mainV03_56_B.Thrust_m;

  /* InitialCondition: '<S4>/IC10' */
  if ((mainV03_56_DW.IC10_FirstOutputTime == (rtMinusInf)) ||
      (mainV03_56_DW.IC10_FirstOutputTime == mainV03_56_M->Timing.t[0])) {
    mainV03_56_DW.IC10_FirstOutputTime = mainV03_56_M->Timing.t[0];
    rtb_sincos_o2_p[0] = mainV03_56_P.IC10_Value[0];
    rtb_sincos_o2_p[1] = mainV03_56_P.IC10_Value[1];
    rtb_sincos_o2_p[2] = mainV03_56_P.IC10_Value[2];
  } else {
    rtb_sincos_o2_p[0] = mainV03_56_B.total_thrust[0];
    rtb_sincos_o2_p[1] = mainV03_56_B.total_thrust[1];
    rtb_sincos_o2_p[2] = mainV03_56_B.total_thrust[2];
  }

  /* End of InitialCondition: '<S4>/IC10' */

  /* Sum: '<S4>/Sum' */
  mainV03_56_B.Sum_ma[0] = ((rtb_Shiftright[0] + mainV03_56_B.envData.Fgravity[0])
    + mainV03_56_B.envData.Fground[0]) + rtb_sincos_o2_p[0];
  mainV03_56_B.Sum_ma[1] = ((rtb_Shiftright[1] + mainV03_56_B.envData.Fgravity[1])
    + mainV03_56_B.envData.Fground[1]) + rtb_sincos_o2_p[1];
  mainV03_56_B.Sum_ma[2] = ((rtb_Shiftright[2] + mainV03_56_B.envData.Fgravity[2])
    + mainV03_56_B.envData.Fground[2]) + rtb_sincos_o2_p[2];

  /* InitialCondition: '<S4>/IC6' */
  if ((mainV03_56_DW.IC6_FirstOutputTime == (rtMinusInf)) ||
      (mainV03_56_DW.IC6_FirstOutputTime == mainV03_56_M->Timing.t[0])) {
    mainV03_56_DW.IC6_FirstOutputTime = mainV03_56_M->Timing.t[0];
    rtb_Shiftright[0] = mainV03_56_P.IC6_Value[0];
    rtb_Shiftright[1] = mainV03_56_P.IC6_Value[1];
    rtb_Shiftright[2] = mainV03_56_P.IC6_Value[2];
  } else {
    rtb_Shiftright[0] = mainV03_56_B.Sum_ma[0];
    rtb_Shiftright[1] = mainV03_56_B.Sum_ma[1];
    rtb_Shiftright[2] = mainV03_56_B.Sum_ma[2];
  }

  /* End of InitialCondition: '<S4>/IC6' */

  /* Product: '<S309>/j x k' */
  mainV03_56_B.jxk_j = mainV03_56_B.ubvbwb[1] * mainV03_56_B.pqr[2];

  /* Product: '<S309>/k x i' */
  mainV03_56_B.kxi_k = mainV03_56_B.ubvbwb[2] * mainV03_56_B.pqr[0];

  /* Product: '<S309>/i x j' */
  mainV03_56_B.ixj_f = mainV03_56_B.ubvbwb[0] * mainV03_56_B.pqr[1];

  /* Product: '<S310>/k x j' */
  mainV03_56_B.kxj_n = mainV03_56_B.ubvbwb[2] * mainV03_56_B.pqr[1];

  /* Product: '<S310>/i x k' */
  mainV03_56_B.ixk_bi = mainV03_56_B.ubvbwb[0] * mainV03_56_B.pqr[2];

  /* Product: '<S310>/j x i' */
  mainV03_56_B.jxi_dx = mainV03_56_B.ubvbwb[1] * mainV03_56_B.pqr[0];

  /* Sum: '<S294>/Sum' */
  mainV03_56_B.Sum_as[0] = mainV03_56_B.jxk_j - mainV03_56_B.kxj_n;
  mainV03_56_B.Sum_as[1] = mainV03_56_B.kxi_k - mainV03_56_B.ixk_bi;
  mainV03_56_B.Sum_as[2] = mainV03_56_B.ixj_f - mainV03_56_B.jxi_dx;

  /* Sum: '<S293>/Sum' */
  mainV03_56_B.Sum_ny[0] = rtb_Shiftright[0] + mainV03_56_B.SumofElements[0];

  /* Product: '<S243>/Product' incorporates:
   *  Constant: '<S246>/mass (kg)'
   */
  mainV03_56_B.Product_km[0] = mainV03_56_B.Sum_ny[0] /
    mainV03_56_P.masskg_Value;

  /* Sum: '<S243>/Sum2' */
  mainV03_56_B.Sum2_af[0] = mainV03_56_B.Product_km[0] + mainV03_56_B.Sum_as[0];

  /* Sum: '<S293>/Sum' */
  mainV03_56_B.Sum_ny[1] = rtb_Shiftright[1] + mainV03_56_B.SumofElements[1];

  /* Product: '<S243>/Product' incorporates:
   *  Constant: '<S246>/mass (kg)'
   */
  mainV03_56_B.Product_km[1] = mainV03_56_B.Sum_ny[1] /
    mainV03_56_P.masskg_Value;

  /* Sum: '<S243>/Sum2' */
  mainV03_56_B.Sum2_af[1] = mainV03_56_B.Product_km[1] + mainV03_56_B.Sum_as[1];

  /* Sum: '<S293>/Sum' */
  mainV03_56_B.Sum_ny[2] = rtb_Shiftright[2] + mainV03_56_B.SumofElements[2];

  /* Product: '<S243>/Product' incorporates:
   *  Constant: '<S246>/mass (kg)'
   */
  mainV03_56_B.Product_km[2] = mainV03_56_B.Sum_ny[2] /
    mainV03_56_P.masskg_Value;

  /* Sum: '<S243>/Sum2' */
  mainV03_56_B.Sum2_af[2] = mainV03_56_B.Product_km[2] + mainV03_56_B.Sum_as[2];

  /* InitialCondition: '<S4>/IC' */
  if ((mainV03_56_DW.IC_FirstOutputTime_n == (rtMinusInf)) ||
      (mainV03_56_DW.IC_FirstOutputTime_n == mainV03_56_M->Timing.t[0])) {
    mainV03_56_DW.IC_FirstOutputTime_n = mainV03_56_M->Timing.t[0];
    mainV03_56_B.dOmega_body_k[0] = mainV03_56_P.IC_Value_l[0];
    mainV03_56_B.dOmega_body_k[1] = mainV03_56_P.IC_Value_l[1];
    mainV03_56_B.dOmega_body_k[2] = mainV03_56_P.IC_Value_l[2];
  } else {
    mainV03_56_B.dOmega_body_k[0] = mainV03_56_B.Product2_ex[0];
    mainV03_56_B.dOmega_body_k[1] = mainV03_56_B.Product2_ex[1];
    mainV03_56_B.dOmega_body_k[2] = mainV03_56_B.Product2_ex[2];
  }

  /* End of InitialCondition: '<S4>/IC' */

  /* InitialCondition: '<S4>/IC1' */
  if ((mainV03_56_DW.IC1_FirstOutputTime_i == (rtMinusInf)) ||
      (mainV03_56_DW.IC1_FirstOutputTime_i == mainV03_56_M->Timing.t[0])) {
    mainV03_56_DW.IC1_FirstOutputTime_i = mainV03_56_M->Timing.t[0];
    mainV03_56_B.Accel_body_a[0] = mainV03_56_P.IC1_Value_n[0];
    mainV03_56_B.Accel_body_a[1] = mainV03_56_P.IC1_Value_n[1];
    mainV03_56_B.Accel_body_a[2] = mainV03_56_P.IC1_Value_n[2];
  } else {
    mainV03_56_B.Accel_body_a[0] = mainV03_56_B.Sum2_af[0];
    mainV03_56_B.Accel_body_a[1] = mainV03_56_B.Sum2_af[1];
    mainV03_56_B.Accel_body_a[2] = mainV03_56_B.Sum2_af[2];
  }

  /* End of InitialCondition: '<S4>/IC1' */

  /* Interpolation_n-D: '<S353>/Current Interpolation Using Prelookup' */
  frac_1d[0] = mainV03_56_B.VoltagePrelookup_o2;
  frac_1d[1] = mainV03_56_B.HeightPrelookup_o2;
  frac_1d[2] = mainV03_56_B.VflightPrelookup_o2;
  bpIndex_1d[0] = mainV03_56_B.VoltagePrelookup_o1;
  bpIndex_1d[1] = mainV03_56_B.HeightPrelookup_o1;
  bpIndex_1d[2] = mainV03_56_B.VflightPrelookup_o1;
  mainV03_56_B.Current = intrp3d_la_pw(bpIndex_1d, frac_1d,
    mainV03_56_P.CurrentInterpolationUsingPrelookup_Table,
    mainV03_56_P.CurrentInterpolationUsingPrelookup_dimSize,
    mainV03_56_P.CurrentInterpolationUsingPrelookup_maxIndex);

  /* Interpolation_n-D: '<S367>/Current Interpolation Using Prelookup' */
  frac_1e[0] = mainV03_56_B.VoltagePrelookup_o2_j;
  frac_1e[1] = mainV03_56_B.HeightPrelookup_o2_k;
  frac_1e[2] = mainV03_56_B.VflightPrelookup_o2_d;
  bpIndex_1e[0] = mainV03_56_B.VoltagePrelookup_o1_a;
  bpIndex_1e[1] = mainV03_56_B.HeightPrelookup_o1_g;
  bpIndex_1e[2] = mainV03_56_B.VflightPrelookup_o1_e;
  mainV03_56_B.Current_d = intrp3d_la_pw(bpIndex_1e, frac_1e,
    mainV03_56_P.CurrentInterpolationUsingPrelookup_Table_l,
    mainV03_56_P.CurrentInterpolationUsingPrelookup_dimSize_e,
    mainV03_56_P.CurrentInterpolationUsingPrelookup_maxIndex_h);

  /* Interpolation_n-D: '<S381>/Current Interpolation Using Prelookup' */
  frac_1f[0] = mainV03_56_B.VoltagePrelookup_o2_a;
  frac_1f[1] = mainV03_56_B.HeightPrelookup_o2_j;
  frac_1f[2] = mainV03_56_B.VflightPrelookup_o2_f;
  bpIndex_1f[0] = mainV03_56_B.VoltagePrelookup_o1_l;
  bpIndex_1f[1] = mainV03_56_B.HeightPrelookup_o1_i;
  bpIndex_1f[2] = mainV03_56_B.VflightPrelookup_o1_m;
  mainV03_56_B.Current_b = intrp3d_la_pw(bpIndex_1f, frac_1f,
    mainV03_56_P.CurrentInterpolationUsingPrelookup_Table_e,
    mainV03_56_P.CurrentInterpolationUsingPrelookup_dimSize_p,
    mainV03_56_P.CurrentInterpolationUsingPrelookup_maxIndex_o);

  /* Interpolation_n-D: '<S395>/Current Interpolation Using Prelookup' */
  frac_1g[0] = mainV03_56_B.VoltagePrelookup_o2_k;
  frac_1g[1] = mainV03_56_B.HeightPrelookup_o2_e;
  frac_1g[2] = mainV03_56_B.VflightPrelookup_o2_fy;
  bpIndex_1g[0] = mainV03_56_B.VoltagePrelookup_o1_i;
  bpIndex_1g[1] = mainV03_56_B.HeightPrelookup_o1_d;
  bpIndex_1g[2] = mainV03_56_B.VflightPrelookup_o1_a;
  mainV03_56_B.Current_m = intrp3d_la_pw(bpIndex_1g, frac_1g,
    mainV03_56_P.CurrentInterpolationUsingPrelookup_Table_k,
    mainV03_56_P.CurrentInterpolationUsingPrelookup_dimSize_k,
    mainV03_56_P.CurrentInterpolationUsingPrelookup_maxIndex_c);

  /* Interpolation_n-D: '<S409>/Current Interpolation Using Prelookup' */
  frac_1h[0] = mainV03_56_B.VoltagePrelookup_o2_f;
  frac_1h[1] = mainV03_56_B.HeightPrelookup_o2_o;
  frac_1h[2] = mainV03_56_B.VflightPrelookup_o2_m;
  bpIndex_1h[0] = mainV03_56_B.VoltagePrelookup_o1_a2;
  bpIndex_1h[1] = mainV03_56_B.HeightPrelookup_o1_l;
  bpIndex_1h[2] = mainV03_56_B.VflightPrelookup_o1_g;
  mainV03_56_B.Current_h = intrp3d_la_pw(bpIndex_1h, frac_1h,
    mainV03_56_P.CurrentInterpolationUsingPrelookup_Table_lx,
    mainV03_56_P.CurrentInterpolationUsingPrelookup_dimSize_i,
    mainV03_56_P.CurrentInterpolationUsingPrelookup_maxIndex_a);

  /* Sum: '<S247>/Sum5' */
  mainV03_56_B.total_Im = (((mainV03_56_B.Current + mainV03_56_B.Current_d) +
    mainV03_56_B.Current_b) + mainV03_56_B.Current_m) + mainV03_56_B.Current_h;
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Sum: '<S464>/Sum' incorporates:
     *  Constant: '<S464>/Constant'
     *  UnitDelay: '<S464>/Unit Delay'
     */
    mainV03_56_B.Sum_db = mainV03_56_P.Constant_Value_is -
      mainV03_56_DW.UnitDelay_DSTATE;
  }

  /* Trigonometry: '<S528>/sincos' */
  rtb_Shiftright[0] = sin(mainV03_56_B.phithetapsi_m[0]);
  rtb_sincos_o2_p[0] = cos(mainV03_56_B.phithetapsi_m[0]);
  rtb_sincos_o2_p[1] = cos(mainV03_56_B.phithetapsi_m[1]);
  rtb_sincos_o2_p[2] = cos(mainV03_56_B.phithetapsi_m[2]);

  /* Fcn: '<S528>/phidot' incorporates:
   *  Trigonometry: '<S528>/sincos'
   */
  mainV03_56_B.phidot_p = (mainV03_56_B.Saturation_a[1] * rtb_Shiftright[0] +
    mainV03_56_B.Saturation_a[2] * rtb_sincos_o2_p[0]) * (sin
    (mainV03_56_B.phithetapsi_m[1]) / rtb_sincos_o2_p[1]) +
    mainV03_56_B.Saturation_a[0];

  /* Fcn: '<S528>/psidot' */
  mainV03_56_B.psidot_k = (mainV03_56_B.Saturation_a[1] * rtb_Shiftright[0] +
    mainV03_56_B.Saturation_a[2] * rtb_sincos_o2_p[0]) / rtb_sincos_o2_p[1];

  /* Fcn: '<S528>/thetadot' */
  mainV03_56_B.thetadot_i = mainV03_56_B.Saturation_a[1] * rtb_sincos_o2_p[0] -
    mainV03_56_B.Saturation_a[2] * rtb_Shiftright[0];
}

/* Model update function */
void mainV03_56_update(void)
{
  real_T *lastU;
  real_T (*lastU_0)[3];

  /* Update for Derivative: '<S4>/Derivative' */
  if (mainV03_56_DW.TimeStampA == (rtInf)) {
    mainV03_56_DW.TimeStampA = mainV03_56_M->Timing.t[0];
    lastU = &mainV03_56_DW.LastUAtTimeA;
  } else if (mainV03_56_DW.TimeStampB == (rtInf)) {
    mainV03_56_DW.TimeStampB = mainV03_56_M->Timing.t[0];
    lastU = &mainV03_56_DW.LastUAtTimeB;
  } else if (mainV03_56_DW.TimeStampA < mainV03_56_DW.TimeStampB) {
    mainV03_56_DW.TimeStampA = mainV03_56_M->Timing.t[0];
    lastU = &mainV03_56_DW.LastUAtTimeA;
  } else {
    mainV03_56_DW.TimeStampB = mainV03_56_M->Timing.t[0];
    lastU = &mainV03_56_DW.LastUAtTimeB;
  }

  *lastU = mainV03_56_B.alpha_j;

  /* End of Update for Derivative: '<S4>/Derivative' */

  /* Update for RateLimiter: '<S4>/Rate Limiter' */
  if (mainV03_56_DW.LastMajorTimeA == (rtInf)) {
    mainV03_56_DW.LastMajorTimeA = mainV03_56_M->Timing.t[0];
    mainV03_56_DW.PrevYA = mainV03_56_B.alpha_dot;
  } else if (mainV03_56_DW.LastMajorTimeB == (rtInf)) {
    mainV03_56_DW.LastMajorTimeB = mainV03_56_M->Timing.t[0];
    mainV03_56_DW.PrevYB = mainV03_56_B.alpha_dot;
  } else if (mainV03_56_DW.LastMajorTimeA < mainV03_56_DW.LastMajorTimeB) {
    mainV03_56_DW.LastMajorTimeA = mainV03_56_M->Timing.t[0];
    mainV03_56_DW.PrevYA = mainV03_56_B.alpha_dot;
  } else {
    mainV03_56_DW.LastMajorTimeB = mainV03_56_M->Timing.t[0];
    mainV03_56_DW.PrevYB = mainV03_56_B.alpha_dot;
  }

  /* End of Update for RateLimiter: '<S4>/Rate Limiter' */

  /* Update for Derivative: '<S4>/Derivative1' */
  if (mainV03_56_DW.TimeStampA_i == (rtInf)) {
    mainV03_56_DW.TimeStampA_i = mainV03_56_M->Timing.t[0];
    lastU = &mainV03_56_DW.LastUAtTimeA_o;
  } else if (mainV03_56_DW.TimeStampB_a == (rtInf)) {
    mainV03_56_DW.TimeStampB_a = mainV03_56_M->Timing.t[0];
    lastU = &mainV03_56_DW.LastUAtTimeB_j;
  } else if (mainV03_56_DW.TimeStampA_i < mainV03_56_DW.TimeStampB_a) {
    mainV03_56_DW.TimeStampA_i = mainV03_56_M->Timing.t[0];
    lastU = &mainV03_56_DW.LastUAtTimeA_o;
  } else {
    mainV03_56_DW.TimeStampB_a = mainV03_56_M->Timing.t[0];
    lastU = &mainV03_56_DW.LastUAtTimeB_j;
  }

  *lastU = mainV03_56_B.beta_l;

  /* End of Update for Derivative: '<S4>/Derivative1' */

  /* Update for RateLimiter: '<S4>/Rate Limiter1' */
  if (mainV03_56_DW.LastMajorTimeA_f == (rtInf)) {
    mainV03_56_DW.LastMajorTimeA_f = mainV03_56_M->Timing.t[0];
    mainV03_56_DW.PrevYA_o = mainV03_56_B.beta_dot;
  } else if (mainV03_56_DW.LastMajorTimeB_d == (rtInf)) {
    mainV03_56_DW.LastMajorTimeB_d = mainV03_56_M->Timing.t[0];
    mainV03_56_DW.PrevYB_d = mainV03_56_B.beta_dot;
  } else if (mainV03_56_DW.LastMajorTimeA_f < mainV03_56_DW.LastMajorTimeB_d) {
    mainV03_56_DW.LastMajorTimeA_f = mainV03_56_M->Timing.t[0];
    mainV03_56_DW.PrevYA_o = mainV03_56_B.beta_dot;
  } else {
    mainV03_56_DW.LastMajorTimeB_d = mainV03_56_M->Timing.t[0];
    mainV03_56_DW.PrevYB_d = mainV03_56_B.beta_dot;
  }

  /* End of Update for RateLimiter: '<S4>/Rate Limiter1' */
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Update for Memory: '<S112>/otime' */
    mainV03_56_DW.otime_PreviousInput = mainV03_56_B.Sum_p;

    /* Update for Memory: '<S111>/olon' */
    mainV03_56_DW.olon_PreviousInput = mainV03_56_B.u80deg;

    /* Update for Memory: '<S110>/olat' */
    mainV03_56_DW.olat_PreviousInput = mainV03_56_B.u0deg;

    /* Update for Memory: '<S110>/oalt' */
    mainV03_56_DW.oalt_PreviousInput = mainV03_56_B.Gain;
  }

  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[2] == 0) {
    /* Update for RandomNumber: '<S65>/White Noise' */
    mainV03_56_DW.NextOutput = rt_nrand_Upu32_Yd_f_pw_snf
      (&mainV03_56_DW.RandSeed) * mainV03_56_P.WhiteNoise_StdDev +
      mainV03_56_P.WhiteNoise_Mean;

    /* Update for RandomNumber: '<S66>/White Noise' */
    mainV03_56_DW.NextOutput_a[0] = rt_nrand_Upu32_Yd_f_pw_snf
      (&mainV03_56_DW.RandSeed_n[0]) * mainV03_56_P.WhiteNoise_StdDev_c +
      mainV03_56_P.WhiteNoise_Mean_o;
    mainV03_56_DW.NextOutput_a[1] = rt_nrand_Upu32_Yd_f_pw_snf
      (&mainV03_56_DW.RandSeed_n[1]) * mainV03_56_P.WhiteNoise_StdDev_c +
      mainV03_56_P.WhiteNoise_Mean_o;
    mainV03_56_DW.NextOutput_a[2] = rt_nrand_Upu32_Yd_f_pw_snf
      (&mainV03_56_DW.RandSeed_n[2]) * mainV03_56_P.WhiteNoise_StdDev_c +
      mainV03_56_P.WhiteNoise_Mean_o;
  }

  /* Update for Enabled SubSystem: '<S55>/Hpgw' incorporates:
   *  Update for EnablePort: '<S67>/Enable'
   */
  if (mainV03_56_DW.Hpgw_MODE && (rtmIsMajorTimeStep(mainV03_56_M) &&
       mainV03_56_M->Timing.TaskCounters.TID[2] == 0)) {
    /* Update for UnitDelay: '<S67>/Unit Delay' */
    mainV03_56_DW.UnitDelay_DSTATE_c[0] = mainV03_56_B.Sum_ax[0];
    mainV03_56_DW.UnitDelay_DSTATE_c[1] = mainV03_56_B.Sum_ax[1];
  }

  /* End of Update for SubSystem: '<S55>/Hpgw' */

  /* Update for Enabled SubSystem: '<S56>/Hwgw(z)' incorporates:
   *  Update for EnablePort: '<S72>/Enable'
   */
  if (mainV03_56_DW.Hwgwz_MODE && (rtmIsMajorTimeStep(mainV03_56_M) &&
       mainV03_56_M->Timing.TaskCounters.TID[2] == 0)) {
    /* Update for UnitDelay: '<S72>/Unit Delay' */
    mainV03_56_DW.UnitDelay_DSTATE_l[0] = mainV03_56_B.Sum_da[0];
    mainV03_56_DW.UnitDelay_DSTATE_l[1] = mainV03_56_B.Sum_da[1];
  }

  /* End of Update for SubSystem: '<S56>/Hwgw(z)' */

  /* Update for Enabled SubSystem: '<S55>/Hqgw' incorporates:
   *  Update for EnablePort: '<S68>/Enable'
   */
  if (mainV03_56_DW.Hqgw_MODE && (rtmIsMajorTimeStep(mainV03_56_M) &&
       mainV03_56_M->Timing.TaskCounters.TID[2] == 0)) {
    /* Update for UnitDelay: '<S68>/Unit Delay' */
    mainV03_56_DW.UnitDelay_DSTATE_fj[0] = mainV03_56_B.Sum1_gx[0];

    /* Update for UnitDelay: '<S68>/Unit Delay1' */
    mainV03_56_DW.UnitDelay1_DSTATE_j[0] = mainV03_56_B.Sum_da[0];

    /* Update for UnitDelay: '<S68>/Unit Delay' */
    mainV03_56_DW.UnitDelay_DSTATE_fj[1] = mainV03_56_B.Sum1_gx[1];

    /* Update for UnitDelay: '<S68>/Unit Delay1' */
    mainV03_56_DW.UnitDelay1_DSTATE_j[1] = mainV03_56_B.Sum_da[1];
  }

  /* End of Update for SubSystem: '<S55>/Hqgw' */

  /* Update for Enabled SubSystem: '<S56>/Hvgw(z)' incorporates:
   *  Update for EnablePort: '<S71>/Enable'
   */
  if (mainV03_56_DW.Hvgwz_MODE && (rtmIsMajorTimeStep(mainV03_56_M) &&
       mainV03_56_M->Timing.TaskCounters.TID[2] == 0)) {
    /* Update for UnitDelay: '<S71>/Unit Delay' */
    mainV03_56_DW.UnitDelay_DSTATE_o[0] = mainV03_56_B.Sum_nt[0];
    mainV03_56_DW.UnitDelay_DSTATE_o[1] = mainV03_56_B.Sum_nt[1];
  }

  /* End of Update for SubSystem: '<S56>/Hvgw(z)' */

  /* Update for Enabled SubSystem: '<S55>/Hrgw' incorporates:
   *  Update for EnablePort: '<S69>/Enable'
   */
  if (mainV03_56_DW.Hrgw_MODE && (rtmIsMajorTimeStep(mainV03_56_M) &&
       mainV03_56_M->Timing.TaskCounters.TID[2] == 0)) {
    /* Update for UnitDelay: '<S69>/Unit Delay' */
    mainV03_56_DW.UnitDelay_DSTATE_ju[0] = mainV03_56_B.Sum1_p2[0];

    /* Update for UnitDelay: '<S69>/Unit Delay1' */
    mainV03_56_DW.UnitDelay1_DSTATE_ml[0] = mainV03_56_B.Sum_nt[0];

    /* Update for UnitDelay: '<S69>/Unit Delay' */
    mainV03_56_DW.UnitDelay_DSTATE_ju[1] = mainV03_56_B.Sum1_p2[1];

    /* Update for UnitDelay: '<S69>/Unit Delay1' */
    mainV03_56_DW.UnitDelay1_DSTATE_ml[1] = mainV03_56_B.Sum_nt[1];
  }

  /* End of Update for SubSystem: '<S55>/Hrgw' */

  /* Update for Enabled SubSystem: '<S56>/Hugw(z)' incorporates:
   *  Update for EnablePort: '<S70>/Enable'
   */
  if (mainV03_56_DW.Hugwz_MODE && (rtmIsMajorTimeStep(mainV03_56_M) &&
       mainV03_56_M->Timing.TaskCounters.TID[2] == 0)) {
    /* Update for UnitDelay: '<S70>/Unit Delay' */
    mainV03_56_DW.UnitDelay_DSTATE_p[0] = mainV03_56_B.Sum_eq[0];
    mainV03_56_DW.UnitDelay_DSTATE_p[1] = mainV03_56_B.Sum_eq[1];
  }

  /* End of Update for SubSystem: '<S56>/Hugw(z)' */

  /* Update for Derivative: '<S498>/Derivative' */
  if (mainV03_56_DW.TimeStampA_p == (rtInf)) {
    mainV03_56_DW.TimeStampA_p = mainV03_56_M->Timing.t[0];
    lastU_0 = (real_T (*)[3])mainV03_56_DW.LastUAtTimeA_g;
  } else if (mainV03_56_DW.TimeStampB_k == (rtInf)) {
    mainV03_56_DW.TimeStampB_k = mainV03_56_M->Timing.t[0];
    lastU_0 = (real_T (*)[3])mainV03_56_DW.LastUAtTimeB_a;
  } else if (mainV03_56_DW.TimeStampA_p < mainV03_56_DW.TimeStampB_k) {
    mainV03_56_DW.TimeStampA_p = mainV03_56_M->Timing.t[0];
    lastU_0 = (real_T (*)[3])mainV03_56_DW.LastUAtTimeA_g;
  } else {
    mainV03_56_DW.TimeStampB_k = mainV03_56_M->Timing.t[0];
    lastU_0 = (real_T (*)[3])mainV03_56_DW.LastUAtTimeB_a;
  }

  (*lastU_0)[0] = mainV03_56_B.X_nedMeas[0];
  (*lastU_0)[1] = mainV03_56_B.X_nedMeas[1];
  (*lastU_0)[2] = mainV03_56_B.X_nedMeas[2];

  /* End of Update for Derivative: '<S498>/Derivative' */
  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[2] == 0) {
    /* Update for RandomNumber: '<S534>/White Noise' */
    mainV03_56_DW.NextOutput_e[0] = rt_nrand_Upu32_Yd_f_pw_snf
      (&mainV03_56_DW.RandSeed_g[0]) * mainV03_56_P.WhiteNoise_StdDev_n +
      mainV03_56_P.WhiteNoise_Mean_d;

    /* Update for RandomNumber: '<S548>/White Noise' */
    mainV03_56_DW.NextOutput_b[0] = rt_nrand_Upu32_Yd_f_pw_snf
      (&mainV03_56_DW.RandSeed_j[0]) * mainV03_56_P.WhiteNoise_StdDev_l +
      mainV03_56_P.WhiteNoise_Mean_n;

    /* Update for RandomNumber: '<S534>/White Noise' */
    mainV03_56_DW.NextOutput_e[1] = rt_nrand_Upu32_Yd_f_pw_snf
      (&mainV03_56_DW.RandSeed_g[1]) * mainV03_56_P.WhiteNoise_StdDev_n +
      mainV03_56_P.WhiteNoise_Mean_d;

    /* Update for RandomNumber: '<S548>/White Noise' */
    mainV03_56_DW.NextOutput_b[1] = rt_nrand_Upu32_Yd_f_pw_snf
      (&mainV03_56_DW.RandSeed_j[1]) * mainV03_56_P.WhiteNoise_StdDev_l +
      mainV03_56_P.WhiteNoise_Mean_n;

    /* Update for RandomNumber: '<S534>/White Noise' */
    mainV03_56_DW.NextOutput_e[2] = rt_nrand_Upu32_Yd_f_pw_snf
      (&mainV03_56_DW.RandSeed_g[2]) * mainV03_56_P.WhiteNoise_StdDev_n +
      mainV03_56_P.WhiteNoise_Mean_d;

    /* Update for RandomNumber: '<S548>/White Noise' */
    mainV03_56_DW.NextOutput_b[2] = rt_nrand_Upu32_Yd_f_pw_snf
      (&mainV03_56_DW.RandSeed_j[2]) * mainV03_56_P.WhiteNoise_StdDev_l +
      mainV03_56_P.WhiteNoise_Mean_n;
  }

  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[1] == 0) {
    /* Update for UnitDelay: '<S464>/Unit Delay' */
    mainV03_56_DW.UnitDelay_DSTATE = mainV03_56_B.Sum_db;
  }

  if (rtmIsMajorTimeStep(mainV03_56_M)) {
    rt_ertODEUpdateContinuousStates(&mainV03_56_M->solverInfo);
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++mainV03_56_M->Timing.clockTick0)) {
    ++mainV03_56_M->Timing.clockTickH0;
  }

  mainV03_56_M->Timing.t[0] = rtsiGetSolverStopTime(&mainV03_56_M->solverInfo);

  {
    /* Update absolute timer for sample time: [0.02s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++mainV03_56_M->Timing.clockTick1)) {
      ++mainV03_56_M->Timing.clockTickH1;
    }

    mainV03_56_M->Timing.t[1] = mainV03_56_M->Timing.clockTick1 *
      mainV03_56_M->Timing.stepSize1 + mainV03_56_M->Timing.clockTickH1 *
      mainV03_56_M->Timing.stepSize1 * 4294967296.0;
  }

  if (rtmIsMajorTimeStep(mainV03_56_M) &&
      mainV03_56_M->Timing.TaskCounters.TID[2] == 0) {
    /* Update absolute timer for sample time: [0.1s, 0.0s] */
    /* The "clockTick2" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick2"
     * and "Timing.stepSize2". Size of "clockTick2" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick2 and the high bits
     * Timing.clockTickH2. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++mainV03_56_M->Timing.clockTick2)) {
      ++mainV03_56_M->Timing.clockTickH2;
    }

    mainV03_56_M->Timing.t[2] = mainV03_56_M->Timing.clockTick2 *
      mainV03_56_M->Timing.stepSize2 + mainV03_56_M->Timing.clockTickH2 *
      mainV03_56_M->Timing.stepSize2 * 4294967296.0;
  }

  rate_scheduler();
}

/* Derivatives for root system: '<Root>' */
void mainV03_56_derivatives(void)
{
  boolean_T lsat;
  boolean_T usat;
  XDot_mainV03_56_T *_rtXdot;
  _rtXdot = ((XDot_mainV03_56_T *) mainV03_56_M->derivs);

  /* Derivatives for Integrator: '<S291>/phi theta psi' */
  _rtXdot->phithetapsi_CSTATE[0] = mainV03_56_B.phidot;
  _rtXdot->phithetapsi_CSTATE[1] = mainV03_56_B.thetadot;
  _rtXdot->phithetapsi_CSTATE[2] = mainV03_56_B.psidot;

  /* Derivatives for Integrator: '<S243>/xe,ye,ze' */
  _rtXdot->xeyeze_CSTATE[0] = mainV03_56_B.Product[0];

  /* Derivatives for Integrator: '<S243>/ub,vb,wb' */
  _rtXdot->ubvbwb_CSTATE[0] = mainV03_56_B.Sum2_af[0];

  /* Derivatives for Integrator: '<S243>/p,q,r ' */
  _rtXdot->pqr_CSTATE[0] = mainV03_56_B.Product2_ex[0];

  /* Derivatives for Integrator: '<S243>/xe,ye,ze' */
  _rtXdot->xeyeze_CSTATE[1] = mainV03_56_B.Product[1];

  /* Derivatives for Integrator: '<S243>/ub,vb,wb' */
  _rtXdot->ubvbwb_CSTATE[1] = mainV03_56_B.Sum2_af[1];

  /* Derivatives for Integrator: '<S243>/p,q,r ' */
  _rtXdot->pqr_CSTATE[1] = mainV03_56_B.Product2_ex[1];

  /* Derivatives for Integrator: '<S243>/xe,ye,ze' */
  _rtXdot->xeyeze_CSTATE[2] = mainV03_56_B.Product[2];

  /* Derivatives for Integrator: '<S243>/ub,vb,wb' */
  _rtXdot->ubvbwb_CSTATE[2] = mainV03_56_B.Sum2_af[2];

  /* Derivatives for Integrator: '<S243>/p,q,r ' */
  _rtXdot->pqr_CSTATE[2] = mainV03_56_B.Product2_ex[2];

  /* Derivatives for TransferFcn: '<S250>/Transfer Fcn' */
  _rtXdot->TransferFcn_CSTATE = 0.0;
  _rtXdot->TransferFcn_CSTATE += mainV03_56_P.TransferFcn_A *
    mainV03_56_X.TransferFcn_CSTATE;
  _rtXdot->TransferFcn_CSTATE += mainV03_56_B.dOmega_body_k[0];

  /* Derivatives for TransferFcn: '<S250>/Transfer Fcn1' */
  _rtXdot->TransferFcn1_CSTATE = 0.0;
  _rtXdot->TransferFcn1_CSTATE += mainV03_56_P.TransferFcn1_A *
    mainV03_56_X.TransferFcn1_CSTATE;
  _rtXdot->TransferFcn1_CSTATE += mainV03_56_B.dOmega_body_k[1];

  /* Derivatives for TransferFcn: '<S250>/Transfer Fcn2' */
  _rtXdot->TransferFcn2_CSTATE = 0.0;
  _rtXdot->TransferFcn2_CSTATE += mainV03_56_P.TransferFcn2_A *
    mainV03_56_X.TransferFcn2_CSTATE;
  _rtXdot->TransferFcn2_CSTATE += mainV03_56_B.dOmega_body_k[2];

  /* Derivatives for TransferFcn: '<S250>/Transfer Fcn3' */
  _rtXdot->TransferFcn3_CSTATE = 0.0;
  _rtXdot->TransferFcn3_CSTATE += mainV03_56_P.TransferFcn3_A *
    mainV03_56_X.TransferFcn3_CSTATE;
  _rtXdot->TransferFcn3_CSTATE += mainV03_56_B.Accel_body_a[0];

  /* Derivatives for TransferFcn: '<S250>/Transfer Fcn4' */
  _rtXdot->TransferFcn4_CSTATE = 0.0;
  _rtXdot->TransferFcn4_CSTATE += mainV03_56_P.TransferFcn4_A *
    mainV03_56_X.TransferFcn4_CSTATE;
  _rtXdot->TransferFcn4_CSTATE += mainV03_56_B.Accel_body_a[1];

  /* Derivatives for TransferFcn: '<S250>/Transfer Fcn5' */
  _rtXdot->TransferFcn5_CSTATE = 0.0;
  _rtXdot->TransferFcn5_CSTATE += mainV03_56_P.TransferFcn5_A *
    mainV03_56_X.TransferFcn5_CSTATE;
  _rtXdot->TransferFcn5_CSTATE += mainV03_56_B.Accel_body_a[2];

  /* Derivatives for TransferFcn: '<S249>/Transfer Fcn1' */
  _rtXdot->TransferFcn1_CSTATE_k = 0.0;
  _rtXdot->TransferFcn1_CSTATE_k += mainV03_56_P.TransferFcn1_A_l *
    mainV03_56_X.TransferFcn1_CSTATE_k;
  _rtXdot->TransferFcn1_CSTATE_k += mainV03_56_B.envData.windVelocity[0];

  /* Derivatives for TransferFcn: '<S249>/Transfer Fcn4' */
  _rtXdot->TransferFcn4_CSTATE_l = 0.0;
  _rtXdot->TransferFcn4_CSTATE_l += mainV03_56_P.TransferFcn4_A_f *
    mainV03_56_X.TransferFcn4_CSTATE_l;
  _rtXdot->TransferFcn4_CSTATE_l += mainV03_56_B.envData.windVelocity[1];

  /* Derivatives for TransferFcn: '<S249>/Transfer Fcn5' */
  _rtXdot->TransferFcn5_CSTATE_f = 0.0;
  _rtXdot->TransferFcn5_CSTATE_f += mainV03_56_P.TransferFcn5_A_m *
    mainV03_56_X.TransferFcn5_CSTATE_f;
  _rtXdot->TransferFcn5_CSTATE_f += mainV03_56_B.envData.windVelocity[2];

  /* Derivatives for Integrator: '<S247>/Integrator' */
  _rtXdot->Integrator_CSTATE = mainV03_56_B.total_Im;

  /* Derivatives for Enabled SubSystem: '<S47>/Distance into gust (x)' */
  if (mainV03_56_DW.Distanceintogustx_MODE) {
    /* Derivatives for Integrator: '<S50>/Distance into Gust (x) (Limited to gust length d)' */
    lsat = (mainV03_56_X.DistanceintoGustxLimitedtogustlengthd_CSTATE_a <=
            mainV03_56_P.DistanceintoGustxLimitedtogustlengthd_LowerSat);
    usat = (mainV03_56_X.DistanceintoGustxLimitedtogustlengthd_CSTATE_a >=
            mainV03_56_P.Distanceintogustx_d_m);
    if (((!lsat) && (!usat)) || (lsat && (mainV03_56_B.Airspeed > 0.0)) || (usat
         && (mainV03_56_B.Airspeed < 0.0))) {
      _rtXdot->DistanceintoGustxLimitedtogustlengthd_CSTATE_a =
        mainV03_56_B.Airspeed;
    } else {
      /* in saturation */
      _rtXdot->DistanceintoGustxLimitedtogustlengthd_CSTATE_a = 0.0;
    }

    /* End of Derivatives for Integrator: '<S50>/Distance into Gust (x) (Limited to gust length d)' */
  } else {
    ((XDot_mainV03_56_T *) mainV03_56_M->derivs)
      ->DistanceintoGustxLimitedtogustlengthd_CSTATE_a = 0.0;
  }

  /* End of Derivatives for SubSystem: '<S47>/Distance into gust (x)' */

  /* Derivatives for Enabled SubSystem: '<S47>/Distance into gust (y)' */
  mainV03_56_Distanceintogusty_Deriv(mainV03_56_B.Airspeed,
    &mainV03_56_DW.Distanceintogusty, (P_Distanceintogusty_mainV03_56_T *)
    &mainV03_56_P.Distanceintogusty, &mainV03_56_X.Distanceintogusty,
    &_rtXdot->Distanceintogusty, mainV03_56_P.Distanceintogusty_d_m);

  /* End of Derivatives for SubSystem: '<S47>/Distance into gust (y)' */

  /* Derivatives for Enabled SubSystem: '<S47>/Distance into gust (z)' */
  mainV03_56_Distanceintogusty_Deriv(mainV03_56_B.Airspeed,
    &mainV03_56_DW.Distanceintogustz, (P_Distanceintogusty_mainV03_56_T *)
    &mainV03_56_P.Distanceintogustz, &mainV03_56_X.Distanceintogustz,
    &_rtXdot->Distanceintogustz, mainV03_56_P.Distanceintogustz_d_m);

  /* End of Derivatives for SubSystem: '<S47>/Distance into gust (z)' */

  /* Derivatives for TransferFcn: '<S538>/Transfer Fcn X' */
  _rtXdot->TransferFcnX_CSTATE[0] = 0.0;
  _rtXdot->TransferFcnX_CSTATE[0] += mainV03_56_P.TransferFcnX_A[0] *
    mainV03_56_X.TransferFcnX_CSTATE[0];
  _rtXdot->TransferFcnX_CSTATE[1] = 0.0;
  _rtXdot->TransferFcnX_CSTATE[0] += mainV03_56_P.TransferFcnX_A[1] *
    mainV03_56_X.TransferFcnX_CSTATE[1];
  _rtXdot->TransferFcnX_CSTATE[1] += mainV03_56_X.TransferFcnX_CSTATE[0];
  _rtXdot->TransferFcnX_CSTATE[0] += mainV03_56_B.Sum4_d[0];

  /* Derivatives for TransferFcn: '<S538>/Transfer Fcn Y' */
  _rtXdot->TransferFcnY_CSTATE[0] = 0.0;
  _rtXdot->TransferFcnY_CSTATE[0] += mainV03_56_P.TransferFcnY_A[0] *
    mainV03_56_X.TransferFcnY_CSTATE[0];
  _rtXdot->TransferFcnY_CSTATE[1] = 0.0;
  _rtXdot->TransferFcnY_CSTATE[0] += mainV03_56_P.TransferFcnY_A[1] *
    mainV03_56_X.TransferFcnY_CSTATE[1];
  _rtXdot->TransferFcnY_CSTATE[1] += mainV03_56_X.TransferFcnY_CSTATE[0];
  _rtXdot->TransferFcnY_CSTATE[0] += mainV03_56_B.Sum4_d[1];

  /* Derivatives for TransferFcn: '<S538>/Transfer Fcn Z' */
  _rtXdot->TransferFcnZ_CSTATE[0] = 0.0;
  _rtXdot->TransferFcnZ_CSTATE[0] += mainV03_56_P.TransferFcnZ_A[0] *
    mainV03_56_X.TransferFcnZ_CSTATE[0];
  _rtXdot->TransferFcnZ_CSTATE[1] = 0.0;
  _rtXdot->TransferFcnZ_CSTATE[0] += mainV03_56_P.TransferFcnZ_A[1] *
    mainV03_56_X.TransferFcnZ_CSTATE[1];
  _rtXdot->TransferFcnZ_CSTATE[1] += mainV03_56_X.TransferFcnZ_CSTATE[0];
  _rtXdot->TransferFcnZ_CSTATE[0] += mainV03_56_B.Sum4_d[2];

  /* Derivatives for TransferFcn: '<S550>/Transfer Fcn X' */
  _rtXdot->TransferFcnX_CSTATE_e[0] = 0.0;
  _rtXdot->TransferFcnX_CSTATE_e[0] += mainV03_56_P.TransferFcnX_A_d[0] *
    mainV03_56_X.TransferFcnX_CSTATE_e[0];
  _rtXdot->TransferFcnX_CSTATE_e[1] = 0.0;
  _rtXdot->TransferFcnX_CSTATE_e[0] += mainV03_56_P.TransferFcnX_A_d[1] *
    mainV03_56_X.TransferFcnX_CSTATE_e[1];
  _rtXdot->TransferFcnX_CSTATE_e[1] += mainV03_56_X.TransferFcnX_CSTATE_e[0];
  _rtXdot->TransferFcnX_CSTATE_e[0] += mainV03_56_B.Sum4_a[0];

  /* Derivatives for TransferFcn: '<S550>/Transfer Fcn Y' */
  _rtXdot->TransferFcnY_CSTATE_j[0] = 0.0;
  _rtXdot->TransferFcnY_CSTATE_j[0] += mainV03_56_P.TransferFcnY_A_h[0] *
    mainV03_56_X.TransferFcnY_CSTATE_j[0];
  _rtXdot->TransferFcnY_CSTATE_j[1] = 0.0;
  _rtXdot->TransferFcnY_CSTATE_j[0] += mainV03_56_P.TransferFcnY_A_h[1] *
    mainV03_56_X.TransferFcnY_CSTATE_j[1];
  _rtXdot->TransferFcnY_CSTATE_j[1] += mainV03_56_X.TransferFcnY_CSTATE_j[0];
  _rtXdot->TransferFcnY_CSTATE_j[0] += mainV03_56_B.Sum4_a[1];

  /* Derivatives for TransferFcn: '<S550>/Transfer Fcn Z' */
  _rtXdot->TransferFcnZ_CSTATE_l[0] = 0.0;
  _rtXdot->TransferFcnZ_CSTATE_l[0] += mainV03_56_P.TransferFcnZ_A_p[0] *
    mainV03_56_X.TransferFcnZ_CSTATE_l[0];
  _rtXdot->TransferFcnZ_CSTATE_l[1] = 0.0;
  _rtXdot->TransferFcnZ_CSTATE_l[0] += mainV03_56_P.TransferFcnZ_A_p[1] *
    mainV03_56_X.TransferFcnZ_CSTATE_l[1];
  _rtXdot->TransferFcnZ_CSTATE_l[1] += mainV03_56_X.TransferFcnZ_CSTATE_l[0];
  _rtXdot->TransferFcnZ_CSTATE_l[0] += mainV03_56_B.Sum4_a[2];

  /* Derivatives for Integrator: '<S525>/phi theta psi' */
  _rtXdot->phithetapsi_CSTATE_e[0] = mainV03_56_B.phidot_p;
  _rtXdot->phithetapsi_CSTATE_e[1] = mainV03_56_B.thetadot_i;
  _rtXdot->phithetapsi_CSTATE_e[2] = mainV03_56_B.psidot_k;

  /* Derivatives for Enabled SubSystem: '<S3>/Quadcopter' */
  if (mainV03_56_DW.Quadcopter_MODE) {
    /* Derivatives for Integrator: '<S179>/Filter' */
    _rtXdot->Filter_CSTATE_ci = mainV03_56_B.FilterCoefficient_ct;

    /* Derivatives for Integrator: '<S179>/Integrator' */
    _rtXdot->Integrator_CSTATE_nu = mainV03_56_B.IntegralGain_mk;

    /* Derivatives for Integrator: '<S184>/Integrator' */
    _rtXdot->Integrator_CSTATE_aw = mainV03_56_B.IntegralGain_hx;

    /* Derivatives for Integrator: '<S184>/Filter' */
    _rtXdot->Filter_CSTATE_aw = mainV03_56_B.FilterCoefficient_dn;

    /* Derivatives for Integrator: '<S180>/Integrator' */
    _rtXdot->Integrator_CSTATE_ij = mainV03_56_B.IntegralGain_nx;

    /* Derivatives for Integrator: '<S180>/Filter' */
    _rtXdot->Filter_CSTATE_jg = mainV03_56_B.FilterCoefficient_im;

    /* Derivatives for Integrator: '<S185>/Integrator' */
    _rtXdot->Integrator_CSTATE_hz = mainV03_56_B.IntegralGain_c;

    /* Derivatives for Integrator: '<S185>/Filter' */
    _rtXdot->Filter_CSTATE_l3 = mainV03_56_B.FilterCoefficient_dj;

    /* Derivatives for Integrator: '<S183>/Integrator' */
    _rtXdot->Integrator_CSTATE_pc = mainV03_56_B.IntegralGain_g4;

    /* Derivatives for Integrator: '<S183>/Filter' */
    _rtXdot->Filter_CSTATE_c0 = mainV03_56_B.FilterCoefficient_dr;

    /* Derivatives for Integrator: '<S191>/Integrator' */
    _rtXdot->Integrator_CSTATE_l2 = mainV03_56_B.IntegralGain_pv;

    /* Derivatives for Integrator: '<S191>/Filter' */
    _rtXdot->Filter_CSTATE_a5 = mainV03_56_B.FilterCoefficient_h;

    /* Derivatives for Integrator: '<S187>/Integrator' */
    _rtXdot->Integrator_CSTATE_am = mainV03_56_B.IntegralGain_kc;

    /* Derivatives for Integrator: '<S187>/Filter' */
    _rtXdot->Filter_CSTATE_mg = mainV03_56_B.FilterCoefficient_dv;

    /* Derivatives for Integrator: '<S192>/Integrator' */
    _rtXdot->Integrator_CSTATE_dkb = mainV03_56_B.IntegralGain_pn;

    /* Derivatives for Integrator: '<S192>/Filter' */
    _rtXdot->Filter_CSTATE_oz = mainV03_56_B.FilterCoefficient_cz;

    /* Derivatives for Integrator: '<S188>/Integrator' */
    _rtXdot->Integrator_CSTATE_i5 = mainV03_56_B.IntegralGain_m1;

    /* Derivatives for Integrator: '<S188>/Filter' */
    _rtXdot->Filter_CSTATE_px = mainV03_56_B.FilterCoefficient_ix;

    /* Derivatives for Integrator: '<S182>/Integrator' */
    _rtXdot->Integrator_CSTATE_ng = mainV03_56_B.IntegralGain_l;

    /* Derivatives for Integrator: '<S182>/Filter' */
    _rtXdot->Filter_CSTATE_d = mainV03_56_B.FilterCoefficient_g;

    /* Derivatives for Integrator: '<S190>/Integrator' */
    _rtXdot->Integrator_CSTATE_cp = mainV03_56_B.IntegralGain_an;

    /* Derivatives for Integrator: '<S190>/Filter' */
    _rtXdot->Filter_CSTATE_gi = mainV03_56_B.FilterCoefficient_dc;

    /* Derivatives for Integrator: '<S186>/Integrator' */
    _rtXdot->Integrator_CSTATE_e1 = mainV03_56_B.IntegralGain_ho;

    /* Derivatives for Integrator: '<S186>/Filter' */
    _rtXdot->Filter_CSTATE_pc = mainV03_56_B.FilterCoefficient_mr;

    /* Derivatives for Integrator: '<S189>/Integrator' */
    _rtXdot->Integrator_CSTATE_gw = mainV03_56_B.IntegralGain_df;

    /* Derivatives for Integrator: '<S189>/Filter' */
    _rtXdot->Filter_CSTATE_py = mainV03_56_B.FilterCoefficient_ko;

    /* Derivatives for Integrator: '<S181>/Filter' */
    _rtXdot->Filter_CSTATE_ib = mainV03_56_B.FilterCoefficient_lz;

    /* Derivatives for Integrator: '<S181>/Integrator' */
    _rtXdot->Integrator_CSTATE_l4 = mainV03_56_B.IntegralGain_ke;
  } else {
    {
      real_T *dx;
      int_T i;
      dx = &(((XDot_mainV03_56_T *) mainV03_56_M->derivs)->Filter_CSTATE_ci);
      for (i=0; i < 28; i++) {
        dx[i] = 0.0;
      }
    }
  }

  /* End of Derivatives for SubSystem: '<S3>/Quadcopter' */

  /* Derivatives for Enabled SubSystem: '<S3>/Quadcopter --> Fixed-Wing' */
  if (mainV03_56_DW.QuadcopterFixedWing_MODE) {
    /* Derivatives for Integrator: '<S201>/Integrator' */
    _rtXdot->Integrator_CSTATE_n1 = mainV03_56_B.IntegralGain_m;

    /* Derivatives for Integrator: '<S201>/Filter' */
    _rtXdot->Filter_CSTATE_nt = mainV03_56_B.FilterCoefficient_j;

    /* Derivatives for Integrator: '<S207>/Integrator' */
    _rtXdot->Integrator_CSTATE_e = mainV03_56_B.IntegralGain_p;

    /* Derivatives for Integrator: '<S207>/Filter' */
    _rtXdot->Filter_CSTATE_m = mainV03_56_B.FilterCoefficient_f;

    /* Derivatives for Integrator: '<S206>/Integrator' */
    _rtXdot->Integrator_CSTATE_il = mainV03_56_B.IntegralGain_ay;

    /* Derivatives for Integrator: '<S206>/Filter' */
    _rtXdot->Filter_CSTATE_i = mainV03_56_B.FilterCoefficient_ln;

    /* Derivatives for Integrator: '<S209>/Integrator' */
    _rtXdot->Integrator_CSTATE_mj = mainV03_56_B.IntegralGain_ox;

    /* Derivatives for Integrator: '<S209>/Filter' */
    _rtXdot->Filter_CSTATE_kd = mainV03_56_B.FilterCoefficient_mb;

    /* Derivatives for Integrator: '<S208>/Integrator' */
    _rtXdot->Integrator_CSTATE_dk = mainV03_56_B.IntegralGain_k;

    /* Derivatives for Integrator: '<S208>/Filter' */
    _rtXdot->Filter_CSTATE_ii = mainV03_56_B.FilterCoefficient_n;

    /* Derivatives for Integrator: '<S211>/Integrator' */
    _rtXdot->Integrator_CSTATE_l = mainV03_56_B.IntegralGain_jv;

    /* Derivatives for Integrator: '<S211>/Filter' */
    _rtXdot->Filter_CSTATE_c = mainV03_56_B.FilterCoefficient_nn;

    /* Derivatives for Integrator: '<S210>/Integrator' */
    _rtXdot->Integrator_CSTATE_p0 = mainV03_56_B.IntegralGain_gz;

    /* Derivatives for Integrator: '<S210>/Filter' */
    _rtXdot->Filter_CSTATE_a = mainV03_56_B.FilterCoefficient_d;

    /* Derivatives for Integrator: '<S200>/Filter' */
    _rtXdot->Filter_CSTATE_mr = mainV03_56_B.FilterCoefficient_mw;

    /* Derivatives for Integrator: '<S200>/Integrator' */
    _rtXdot->Integrator_CSTATE_mx = mainV03_56_B.IntegralGain_h;

    /* Derivatives for Integrator: '<S212>/Integrator' */
    _rtXdot->Integrator_CSTATE_g = mainV03_56_B.IntegralGain_bz;

    /* Derivatives for Integrator: '<S212>/Filter' */
    _rtXdot->Filter_CSTATE_o = mainV03_56_B.FilterCoefficient_p;

    /* Derivatives for Integrator: '<S202>/Integrator' */
    _rtXdot->Integrator_CSTATE_o = mainV03_56_B.IntegralGain_e;

    /* Derivatives for Integrator: '<S202>/Filter' */
    _rtXdot->Filter_CSTATE_mk = mainV03_56_B.FilterCoefficient_b;

    /* Derivatives for Integrator: '<S213>/Integrator' */
    _rtXdot->Integrator_CSTATE_io = mainV03_56_B.IntegralGain_n;

    /* Derivatives for Integrator: '<S213>/Filter' */
    _rtXdot->Filter_CSTATE_p = mainV03_56_B.FilterCoefficient_kl;

    /* Derivatives for Integrator: '<S205>/Integrator' */
    _rtXdot->Integrator_CSTATE_m3 = mainV03_56_B.IntegralGain_g;

    /* Derivatives for Integrator: '<S205>/Filter' */
    _rtXdot->Filter_CSTATE_kdq = mainV03_56_B.FilterCoefficient_i;

    /* Derivatives for Integrator: '<S219>/Integrator' */
    _rtXdot->Integrator_CSTATE_gq = mainV03_56_B.IntegralGain_eu;

    /* Derivatives for Integrator: '<S219>/Filter' */
    _rtXdot->Filter_CSTATE_az = mainV03_56_B.FilterCoefficient_mf;

    /* Derivatives for Integrator: '<S215>/Integrator' */
    _rtXdot->Integrator_CSTATE_nb = mainV03_56_B.IntegralGain_ix;

    /* Derivatives for Integrator: '<S215>/Filter' */
    _rtXdot->Filter_CSTATE_f = mainV03_56_B.FilterCoefficient_ps;

    /* Derivatives for Integrator: '<S220>/Integrator' */
    _rtXdot->Integrator_CSTATE_j = mainV03_56_B.IntegralGain_ao;

    /* Derivatives for Integrator: '<S220>/Filter' */
    _rtXdot->Filter_CSTATE_ok = mainV03_56_B.FilterCoefficient_ip;

    /* Derivatives for Integrator: '<S216>/Integrator' */
    _rtXdot->Integrator_CSTATE_h = mainV03_56_B.IntegralGain_kz;

    /* Derivatives for Integrator: '<S216>/Filter' */
    _rtXdot->Filter_CSTATE_or = mainV03_56_B.FilterCoefficient_pt;

    /* Derivatives for Integrator: '<S204>/Integrator' */
    _rtXdot->Integrator_CSTATE_b = mainV03_56_B.IntegralGain_o;

    /* Derivatives for Integrator: '<S204>/Filter' */
    _rtXdot->Filter_CSTATE_j = mainV03_56_B.FilterCoefficient_o5;

    /* Derivatives for Integrator: '<S218>/Integrator' */
    _rtXdot->Integrator_CSTATE_e5 = mainV03_56_B.IntegralGain_da;

    /* Derivatives for Integrator: '<S218>/Filter' */
    _rtXdot->Filter_CSTATE_mkr = mainV03_56_B.FilterCoefficient_m4;

    /* Derivatives for Integrator: '<S214>/Integrator' */
    _rtXdot->Integrator_CSTATE_ni = mainV03_56_B.IntegralGain_d1;

    /* Derivatives for Integrator: '<S214>/Filter' */
    _rtXdot->Filter_CSTATE_l = mainV03_56_B.FilterCoefficient_dq;

    /* Derivatives for Integrator: '<S217>/Integrator' */
    _rtXdot->Integrator_CSTATE_c = mainV03_56_B.IntegralGain_b1;

    /* Derivatives for Integrator: '<S217>/Filter' */
    _rtXdot->Filter_CSTATE_j2 = mainV03_56_B.FilterCoefficient_fx;

    /* Derivatives for Integrator: '<S203>/Filter' */
    _rtXdot->Filter_CSTATE_br = mainV03_56_B.FilterCoefficient_lp;

    /* Derivatives for Integrator: '<S203>/Integrator' */
    _rtXdot->Integrator_CSTATE_a = mainV03_56_B.IntegralGain_fs;
  } else {
    {
      real_T *dx;
      int_T i;
      dx = &(((XDot_mainV03_56_T *) mainV03_56_M->derivs)->Integrator_CSTATE_n1);
      for (i=0; i < 42; i++) {
        dx[i] = 0.0;
      }
    }
  }

  /* End of Derivatives for SubSystem: '<S3>/Quadcopter --> Fixed-Wing' */

  /* Derivatives for Enabled SubSystem: '<S3>/Fixed-Wing Climb' */
  if (mainV03_56_DW.FixedWingClimb_MODE) {
    /* Derivatives for Integrator: '<S176>/Integrator' */
    _rtXdot->Integrator_CSTATE_dg = mainV03_56_B.IntegralGain_cv;

    /* Derivatives for Integrator: '<S176>/Filter' */
    _rtXdot->Filter_CSTATE_k1 = mainV03_56_B.FilterCoefficient_ge;

    /* Derivatives for Integrator: '<S169>/Integrator' */
    _rtXdot->Integrator_CSTATE_ke = mainV03_56_B.IntegralGain_m3;

    /* Derivatives for Integrator: '<S169>/Filter' */
    _rtXdot->Filter_CSTATE_ba = mainV03_56_B.FilterCoefficient_f0;

    /* Derivatives for Integrator: '<S171>/Integrator' */
    _rtXdot->Integrator_CSTATE_ih = mainV03_56_B.IntegralGain_ah;

    /* Derivatives for Integrator: '<S171>/Filter' */
    _rtXdot->Filter_CSTATE_lp = mainV03_56_B.FilterCoefficient_n0;

    /* Derivatives for Integrator: '<S170>/Integrator' */
    _rtXdot->Integrator_CSTATE_nm = mainV03_56_B.IntegralGain_jo;

    /* Derivatives for Integrator: '<S170>/Filter' */
    _rtXdot->Filter_CSTATE_dt = mainV03_56_B.FilterCoefficient_jh;

    /* Derivatives for Integrator: '<S173>/Integrator' */
    _rtXdot->Integrator_CSTATE_og = mainV03_56_B.IntegralGain_pd;

    /* Derivatives for Integrator: '<S173>/Filter' */
    _rtXdot->Filter_CSTATE_iv = mainV03_56_B.FilterCoefficient_d3;

    /* Derivatives for Integrator: '<S172>/Integrator' */
    _rtXdot->Integrator_CSTATE_ht = mainV03_56_B.IntegralGain_fk;

    /* Derivatives for Integrator: '<S172>/Filter' */
    _rtXdot->Filter_CSTATE_k5 = mainV03_56_B.FilterCoefficient_hr;

    /* Derivatives for Integrator: '<S175>/Integrator' */
    _rtXdot->Integrator_CSTATE_cw = mainV03_56_B.IntegralGain_a3;

    /* Derivatives for Integrator: '<S175>/Filter' */
    _rtXdot->Filter_CSTATE_lo = mainV03_56_B.FilterCoefficient_ok;

    /* Derivatives for Integrator: '<S174>/Integrator' */
    _rtXdot->Integrator_CSTATE_ki = mainV03_56_B.IntegralGain_hb;

    /* Derivatives for Integrator: '<S174>/Filter' */
    _rtXdot->Filter_CSTATE_kc = mainV03_56_B.FilterCoefficient_cv;
  } else {
    {
      real_T *dx;
      int_T i;
      dx = &(((XDot_mainV03_56_T *) mainV03_56_M->derivs)->Integrator_CSTATE_dg);
      for (i=0; i < 16; i++) {
        dx[i] = 0.0;
      }
    }
  }

  /* End of Derivatives for SubSystem: '<S3>/Fixed-Wing Climb' */

  /* Derivatives for Enabled SubSystem: '<S3>/TakeOff' */
  if (mainV03_56_DW.TakeOff_MODE) {
    /* Derivatives for Integrator: '<S233>/Integrator' */
    _rtXdot->Integrator_CSTATE_p = mainV03_56_B.IntegralGain;

    /* Derivatives for Integrator: '<S233>/Filter' */
    _rtXdot->Filter_CSTATE = mainV03_56_B.FilterCoefficient;

    /* Derivatives for Integrator: '<S235>/Integrator' */
    _rtXdot->Integrator_CSTATE_i = mainV03_56_B.IntegralGain_i;

    /* Derivatives for Integrator: '<S235>/Filter' */
    _rtXdot->Filter_CSTATE_e = mainV03_56_B.FilterCoefficient_l;

    /* Derivatives for Integrator: '<S234>/Integrator' */
    _rtXdot->Integrator_CSTATE_n = mainV03_56_B.IntegralGain_d;

    /* Derivatives for Integrator: '<S234>/Filter' */
    _rtXdot->Filter_CSTATE_b = mainV03_56_B.FilterCoefficient_c;

    /* Derivatives for Integrator: '<S237>/Integrator' */
    _rtXdot->Integrator_CSTATE_d = mainV03_56_B.IntegralGain_a;

    /* Derivatives for Integrator: '<S237>/Filter' */
    _rtXdot->Filter_CSTATE_g = mainV03_56_B.FilterCoefficient_k;

    /* Derivatives for Integrator: '<S236>/Integrator' */
    _rtXdot->Integrator_CSTATE_k = mainV03_56_B.IntegralGain_f;

    /* Derivatives for Integrator: '<S236>/Filter' */
    _rtXdot->Filter_CSTATE_gp = mainV03_56_B.FilterCoefficient_o;

    /* Derivatives for Integrator: '<S239>/Integrator' */
    _rtXdot->Integrator_CSTATE_m = mainV03_56_B.IntegralGain_b;

    /* Derivatives for Integrator: '<S239>/Filter' */
    _rtXdot->Filter_CSTATE_k = mainV03_56_B.FilterCoefficient_km;

    /* Derivatives for Integrator: '<S238>/Integrator' */
    _rtXdot->Integrator_CSTATE_py = mainV03_56_B.IntegralGain_j;

    /* Derivatives for Integrator: '<S238>/Filter' */
    _rtXdot->Filter_CSTATE_n = mainV03_56_B.FilterCoefficient_m;
  } else {
    {
      real_T *dx;
      int_T i;
      dx = &(((XDot_mainV03_56_T *) mainV03_56_M->derivs)->Integrator_CSTATE_p);
      for (i=0; i < 14; i++) {
        dx[i] = 0.0;
      }
    }
  }

  /* End of Derivatives for SubSystem: '<S3>/TakeOff' */
}

/* Model initialize function */
void mainV03_56_initialize(void)
{
  /* Start for InitialCondition: '<S4>/IC4' */
  mainV03_56_DW.IC4_FirstOutputTime = (rtMinusInf);

  /* Start for InitialCondition: '<S4>/IC11' */
  mainV03_56_DW.IC11_FirstOutputTime = (rtMinusInf);

  /* Start for InitialCondition: '<S4>/IC12' */
  mainV03_56_DW.IC12_FirstOutputTime = (rtMinusInf);

  /* Start for InitialCondition: '<S4>/IC2' */
  mainV03_56_DW.IC2_FirstOutputTime = (rtMinusInf);

  /* Start for InitialCondition: '<S4>/IC5' */
  mainV03_56_DW.IC5_FirstOutputTime = (rtMinusInf);

  /* Start for InitialCondition: '<S4>/IC3' */
  mainV03_56_DW.IC3_FirstOutputTime = (rtMinusInf);

  /* Start for InitialCondition: '<S4>/IC7' */
  mainV03_56_DW.IC7_FirstOutputTime = (rtMinusInf);

  /* Start for InitialCondition: '<S4>/IC8' */
  mainV03_56_DW.IC8_FirstOutputTime = (rtMinusInf);

  /* Start for If: '<S553>/If' */
  mainV03_56_DW.If_ActiveSubsystem = -1;

  /* Start for S-Function (saeroatmos): '<S12>/S-Function' */
  {
    real_T *temp_table = (real_T *) &mainV03_56_DW.SFunction_temp_table[0];
    real_T *pres_table = (real_T *) &mainV03_56_DW.SFunction_pres_table[0];

    /* COESA */
    /*
     * Initialize COESA pressure and temperature tables.
     */
    InitCalcAtmosCOESA( temp_table, pres_table );
  }

  /* Start for InitialCondition: '<S13>/IC' */
  mainV03_56_DW.IC_FirstOutputTime = (rtMinusInf);

  /* Start for InitialCondition: '<S13>/IC1' */
  mainV03_56_DW.IC1_FirstOutputTime = (rtMinusInf);

  /* Start for InitialCondition: '<S13>/IC2' */
  mainV03_56_DW.IC2_FirstOutputTime_d = (rtMinusInf);

  /* Start for InitialCondition: '<S13>/IC3' */
  mainV03_56_DW.IC3_FirstOutputTime_o = (rtMinusInf);

  /* Start for InitialCondition: '<S13>/IC4' */
  mainV03_56_DW.IC4_FirstOutputTime_m = (rtMinusInf);

  /* Start for InitialCondition: '<S13>/IC5' */
  mainV03_56_DW.IC5_FirstOutputTime_d = (rtMinusInf);

  /* Start for If: '<S60>/if Height < Max low altitude  elseif Height > Min isotropic altitude ' */
  mainV03_56_DW.ifHeightMaxlowaltitudeelseifHeightMinisotropicaltitude_ActiveSubsystem
    = -1;

  /* Start for If: '<S61>/if Height < Max low altitude  elseif Height > Min isotropic altitude ' */
  mainV03_56_DW.ifHeightMaxlowaltitudeelseifHeightMinisotropicaltitude_ActiveSubsystem_d
    = -1;

  /* Start for FromWorkspace: '<S463>/FromWs' */
  {
    static real_T pTimeValues0[] = { 0.0, 15.0, 15.0, 20.0, 20.0, 30.0, 30.0,
      35.0, 35.0, 45.0, 45.0, 50.0, 50.0, 500.0 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 2.2, 2.2, 0.0, 0.0, 3.0, 3.0, 0.0, 0.0
    } ;

    mainV03_56_DW.FromWs_PWORK.TimePtr = (void *) pTimeValues0;
    mainV03_56_DW.FromWs_PWORK.DataPtr = (void *) pDataValues0;
    mainV03_56_DW.FromWs_IWORK.PrevIndex = 0;
  }

  /* Level2 S-Function Block: '<S5>/Joystick Input' (joyinput) */
  {
    SimStruct *rts = mainV03_56_M->childSfunctions[0];
    sfcnStart(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Start for Enabled SubSystem: '<S3>/Fixed-Wing - Cruise' */
  /* Start for InitialCondition: '<S163>/IC' */
  mainV03_56_B.Throttle1 = mainV03_56_P.IC_Value;
  mainV03_56_DW.IC_FirstOutputTime_m = true;

  /* Start for InitialCondition: '<S163>/IC1' */
  mainV03_56_B.Throttle2_g = mainV03_56_P.IC1_Value;
  mainV03_56_DW.IC1_FirstOutputTime_c = true;

  /* Start for InitialCondition: '<S163>/IC2' */
  mainV03_56_B.Throttle3_jz = mainV03_56_P.IC2_Value;
  mainV03_56_DW.IC2_FirstOutputTime_do = true;

  /* Start for InitialCondition: '<S163>/IC3' */
  mainV03_56_B.Throttle4_j = mainV03_56_P.IC3_Value;
  mainV03_56_DW.IC3_FirstOutputTime_m = true;

  /* Start for InitialCondition: '<S163>/IC4' */
  mainV03_56_B.Throttle5_a = mainV03_56_P.IC4_Value;
  mainV03_56_DW.IC4_FirstOutputTime_c = true;

  /* End of Start for SubSystem: '<S3>/Fixed-Wing - Cruise' */
  /* Start for InitialCondition: '<S242>/IC1' */
  mainV03_56_DW.IC1_FirstOutputTime_k = (rtMinusInf);

  /* Start for InitialCondition: '<S242>/IC2' */
  mainV03_56_DW.IC2_FirstOutputTime_g = (rtMinusInf);

  /* Start for InitialCondition: '<S4>/IC9' */
  mainV03_56_DW.IC9_FirstOutputTime = (rtMinusInf);

  /* Start for InitialCondition: '<S4>/IC13' */
  mainV03_56_DW.IC13_FirstOutputTime = (rtMinusInf);

  /* Start for InitialCondition: '<S4>/IC10' */
  mainV03_56_DW.IC10_FirstOutputTime = (rtMinusInf);

  /* Start for InitialCondition: '<S4>/IC6' */
  mainV03_56_DW.IC6_FirstOutputTime = (rtMinusInf);

  /* Start for InitialCondition: '<S4>/IC' */
  mainV03_56_DW.IC_FirstOutputTime_n = (rtMinusInf);

  /* Start for InitialCondition: '<S4>/IC1' */
  mainV03_56_DW.IC1_FirstOutputTime_i = (rtMinusInf);
  mainV03_56_PrevZCX.JKFlipFlop_jr.JKFlipFlop_Trig_ZCE = UNINITIALIZED_ZCSIG;
  mainV03_56_PrevZCX.JKFlipFlop_j.JKFlipFlop_Trig_ZCE = UNINITIALIZED_ZCSIG;
  mainV03_56_PrevZCX.JKFlipFlop_f.JKFlipFlop_Trig_ZCE = UNINITIALIZED_ZCSIG;
  mainV03_56_PrevZCX.JKFlipFlop_h.JKFlipFlop_Trig_ZCE = UNINITIALIZED_ZCSIG;
  mainV03_56_PrevZCX.JKFlipFlop_ej.JKFlipFlop_Trig_ZCE = UNINITIALIZED_ZCSIG;
  mainV03_56_PrevZCX.JKFlipFlop_n.JKFlipFlop_Trig_ZCE = UNINITIALIZED_ZCSIG;
  mainV03_56_PrevZCX.JKFlipFlop_b.JKFlipFlop_Trig_ZCE = UNINITIALIZED_ZCSIG;
  mainV03_56_PrevZCX.JKFlipFlop_p1.JKFlipFlop_Trig_ZCE = UNINITIALIZED_ZCSIG;
  mainV03_56_PrevZCX.JKFlipFlop_a.JKFlipFlop_Trig_ZCE = UNINITIALIZED_ZCSIG;
  mainV03_56_PrevZCX.JKFlipFlop_p.JKFlipFlop_Trig_ZCE = UNINITIALIZED_ZCSIG;
  mainV03_56_PrevZCX.JKFlipFlop_e.JKFlipFlop_Trig_ZCE = UNINITIALIZED_ZCSIG;
  mainV03_56_PrevZCX.JKFlipFlop.JKFlipFlop_Trig_ZCE = UNINITIALIZED_ZCSIG;

  {
    uint32_T tseed;
    int32_T r;
    int32_T t;
    real_T y1;

    /* InitializeConditions for Integrator: '<S243>/xe,ye,ze' */
    mainV03_56_X.xeyeze_CSTATE[0] =
      mainV03_56_P.CustomVariableMass6DOFEulerAnglespropio_xme_0[0];

    /* InitializeConditions for Integrator: '<S291>/phi theta psi' */
    mainV03_56_X.phithetapsi_CSTATE[0] =
      mainV03_56_P.CustomVariableMass6DOFEulerAnglespropio_eul_0[0];

    /* InitializeConditions for Integrator: '<S243>/ub,vb,wb' */
    mainV03_56_X.ubvbwb_CSTATE[0] =
      mainV03_56_P.CustomVariableMass6DOFEulerAnglespropio_Vm_0[0];

    /* InitializeConditions for Integrator: '<S243>/p,q,r ' */
    mainV03_56_X.pqr_CSTATE[0] =
      mainV03_56_P.CustomVariableMass6DOFEulerAnglespropio_pm_0[0];

    /* InitializeConditions for Integrator: '<S243>/xe,ye,ze' */
    mainV03_56_X.xeyeze_CSTATE[1] =
      mainV03_56_P.CustomVariableMass6DOFEulerAnglespropio_xme_0[1];

    /* InitializeConditions for Integrator: '<S291>/phi theta psi' */
    mainV03_56_X.phithetapsi_CSTATE[1] =
      mainV03_56_P.CustomVariableMass6DOFEulerAnglespropio_eul_0[1];

    /* InitializeConditions for Integrator: '<S243>/ub,vb,wb' */
    mainV03_56_X.ubvbwb_CSTATE[1] =
      mainV03_56_P.CustomVariableMass6DOFEulerAnglespropio_Vm_0[1];

    /* InitializeConditions for Integrator: '<S243>/p,q,r ' */
    mainV03_56_X.pqr_CSTATE[1] =
      mainV03_56_P.CustomVariableMass6DOFEulerAnglespropio_pm_0[1];

    /* InitializeConditions for Integrator: '<S243>/xe,ye,ze' */
    mainV03_56_X.xeyeze_CSTATE[2] =
      mainV03_56_P.CustomVariableMass6DOFEulerAnglespropio_xme_0[2];

    /* InitializeConditions for Integrator: '<S291>/phi theta psi' */
    mainV03_56_X.phithetapsi_CSTATE[2] =
      mainV03_56_P.CustomVariableMass6DOFEulerAnglespropio_eul_0[2];

    /* InitializeConditions for Integrator: '<S243>/ub,vb,wb' */
    mainV03_56_X.ubvbwb_CSTATE[2] =
      mainV03_56_P.CustomVariableMass6DOFEulerAnglespropio_Vm_0[2];

    /* InitializeConditions for Integrator: '<S243>/p,q,r ' */
    mainV03_56_X.pqr_CSTATE[2] =
      mainV03_56_P.CustomVariableMass6DOFEulerAnglespropio_pm_0[2];

    /* InitializeConditions for TransferFcn: '<S250>/Transfer Fcn' */
    mainV03_56_X.TransferFcn_CSTATE = 0.0;

    /* InitializeConditions for TransferFcn: '<S250>/Transfer Fcn1' */
    mainV03_56_X.TransferFcn1_CSTATE = 0.0;

    /* InitializeConditions for TransferFcn: '<S250>/Transfer Fcn2' */
    mainV03_56_X.TransferFcn2_CSTATE = 0.0;

    /* InitializeConditions for TransferFcn: '<S250>/Transfer Fcn3' */
    mainV03_56_X.TransferFcn3_CSTATE = 0.0;

    /* InitializeConditions for TransferFcn: '<S250>/Transfer Fcn4' */
    mainV03_56_X.TransferFcn4_CSTATE = 0.0;

    /* InitializeConditions for TransferFcn: '<S250>/Transfer Fcn5' */
    mainV03_56_X.TransferFcn5_CSTATE = 0.0;

    /* InitializeConditions for TransferFcn: '<S249>/Transfer Fcn1' */
    mainV03_56_X.TransferFcn1_CSTATE_k = 0.0;

    /* InitializeConditions for TransferFcn: '<S249>/Transfer Fcn4' */
    mainV03_56_X.TransferFcn4_CSTATE_l = 0.0;

    /* InitializeConditions for TransferFcn: '<S249>/Transfer Fcn5' */
    mainV03_56_X.TransferFcn5_CSTATE_f = 0.0;

    /* InitializeConditions for Derivative: '<S4>/Derivative' */
    mainV03_56_DW.TimeStampA = (rtInf);
    mainV03_56_DW.TimeStampB = (rtInf);

    /* InitializeConditions for RateLimiter: '<S4>/Rate Limiter' */
    mainV03_56_DW.LastMajorTimeA = (rtInf);
    mainV03_56_DW.LastMajorTimeB = (rtInf);

    /* InitializeConditions for Derivative: '<S4>/Derivative1' */
    mainV03_56_DW.TimeStampA_i = (rtInf);
    mainV03_56_DW.TimeStampB_a = (rtInf);

    /* InitializeConditions for RateLimiter: '<S4>/Rate Limiter1' */
    mainV03_56_DW.LastMajorTimeA_f = (rtInf);
    mainV03_56_DW.LastMajorTimeB_d = (rtInf);

    /* InitializeConditions for Integrator: '<S247>/Integrator' */
    mainV03_56_X.Integrator_CSTATE = mainV03_56_P.Integrator_IC_evh;

    /* InitializeConditions for Memory: '<S112>/otime' */
    mainV03_56_DW.otime_PreviousInput = mainV03_56_P.otime_X0;

    /* InitializeConditions for Memory: '<S111>/olon' */
    mainV03_56_DW.olon_PreviousInput = mainV03_56_P.olon_X0;

    /* InitializeConditions for Memory: '<S110>/olat' */
    mainV03_56_DW.olat_PreviousInput = mainV03_56_P.olat_X0;

    /* InitializeConditions for Memory: '<S110>/oalt' */
    mainV03_56_DW.oalt_PreviousInput = mainV03_56_P.oalt_X0;

    /* InitializeConditions for RandomNumber: '<S65>/White Noise' */
    y1 = floor(mainV03_56_P.WhiteNoise_Seed);
    if (rtIsNaN(y1) || rtIsInf(y1)) {
      y1 = 0.0;
    } else {
      y1 = fmod(y1, 4.294967296E+9);
    }

    tseed = y1 < 0.0 ? (uint32_T)-(int32_T)(uint32_T)-y1 : (uint32_T)y1;
    r = (int32_T)(tseed >> 16U);
    t = (int32_T)(tseed & 32768U);
    tseed = ((((tseed - ((uint32_T)r << 16U)) + t) << 16U) + t) + r;
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else {
      if (tseed > 2147483646U) {
        tseed = 2147483646U;
      }
    }

    mainV03_56_DW.RandSeed = tseed;
    mainV03_56_DW.NextOutput = rt_nrand_Upu32_Yd_f_pw_snf
      (&mainV03_56_DW.RandSeed) * mainV03_56_P.WhiteNoise_StdDev +
      mainV03_56_P.WhiteNoise_Mean;

    /* End of InitializeConditions for RandomNumber: '<S65>/White Noise' */

    /* InitializeConditions for RandomNumber: '<S66>/White Noise' */
    y1 = floor(mainV03_56_P.WhiteNoise_Seed_f[0]);
    if (rtIsNaN(y1) || rtIsInf(y1)) {
      y1 = 0.0;
    } else {
      y1 = fmod(y1, 4.294967296E+9);
    }

    tseed = y1 < 0.0 ? (uint32_T)-(int32_T)(uint32_T)-y1 : (uint32_T)y1;
    r = (int32_T)(tseed >> 16U);
    t = (int32_T)(tseed & 32768U);
    tseed = ((((tseed - ((uint32_T)r << 16U)) + t) << 16U) + t) + r;
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else {
      if (tseed > 2147483646U) {
        tseed = 2147483646U;
      }
    }

    y1 = rt_nrand_Upu32_Yd_f_pw_snf(&tseed) * mainV03_56_P.WhiteNoise_StdDev_c +
      mainV03_56_P.WhiteNoise_Mean_o;
    mainV03_56_DW.NextOutput_a[0] = y1;
    mainV03_56_DW.RandSeed_n[0] = tseed;
    y1 = floor(mainV03_56_P.WhiteNoise_Seed_f[1]);
    if (rtIsNaN(y1) || rtIsInf(y1)) {
      y1 = 0.0;
    } else {
      y1 = fmod(y1, 4.294967296E+9);
    }

    tseed = y1 < 0.0 ? (uint32_T)-(int32_T)(uint32_T)-y1 : (uint32_T)y1;
    r = (int32_T)(tseed >> 16U);
    t = (int32_T)(tseed & 32768U);
    tseed = ((((tseed - ((uint32_T)r << 16U)) + t) << 16U) + t) + r;
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else {
      if (tseed > 2147483646U) {
        tseed = 2147483646U;
      }
    }

    y1 = rt_nrand_Upu32_Yd_f_pw_snf(&tseed) * mainV03_56_P.WhiteNoise_StdDev_c +
      mainV03_56_P.WhiteNoise_Mean_o;
    mainV03_56_DW.NextOutput_a[1] = y1;
    mainV03_56_DW.RandSeed_n[1] = tseed;
    y1 = floor(mainV03_56_P.WhiteNoise_Seed_f[2]);
    if (rtIsNaN(y1) || rtIsInf(y1)) {
      y1 = 0.0;
    } else {
      y1 = fmod(y1, 4.294967296E+9);
    }

    tseed = y1 < 0.0 ? (uint32_T)-(int32_T)(uint32_T)-y1 : (uint32_T)y1;
    r = (int32_T)(tseed >> 16U);
    t = (int32_T)(tseed & 32768U);
    tseed = ((((tseed - ((uint32_T)r << 16U)) + t) << 16U) + t) + r;
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else {
      if (tseed > 2147483646U) {
        tseed = 2147483646U;
      }
    }

    y1 = rt_nrand_Upu32_Yd_f_pw_snf(&tseed) * mainV03_56_P.WhiteNoise_StdDev_c +
      mainV03_56_P.WhiteNoise_Mean_o;
    mainV03_56_DW.NextOutput_a[2] = y1;
    mainV03_56_DW.RandSeed_n[2] = tseed;

    /* End of InitializeConditions for RandomNumber: '<S66>/White Noise' */

    /* InitializeConditions for Derivative: '<S498>/Derivative' */
    mainV03_56_DW.TimeStampA_p = (rtInf);
    mainV03_56_DW.TimeStampB_k = (rtInf);

    /* InitializeConditions for TransferFcn: '<S538>/Transfer Fcn X' */
    mainV03_56_X.TransferFcnX_CSTATE[0] = 0.0;

    /* InitializeConditions for TransferFcn: '<S538>/Transfer Fcn Y' */
    mainV03_56_X.TransferFcnY_CSTATE[0] = 0.0;

    /* InitializeConditions for TransferFcn: '<S538>/Transfer Fcn Z' */
    mainV03_56_X.TransferFcnZ_CSTATE[0] = 0.0;

    /* InitializeConditions for TransferFcn: '<S538>/Transfer Fcn X' */
    mainV03_56_X.TransferFcnX_CSTATE[1] = 0.0;

    /* InitializeConditions for TransferFcn: '<S538>/Transfer Fcn Y' */
    mainV03_56_X.TransferFcnY_CSTATE[1] = 0.0;

    /* InitializeConditions for TransferFcn: '<S538>/Transfer Fcn Z' */
    mainV03_56_X.TransferFcnZ_CSTATE[1] = 0.0;

    /* InitializeConditions for RandomNumber: '<S534>/White Noise' */
    y1 = floor(mainV03_56_P.WhiteNoise_Seed_n[0]);
    if (rtIsNaN(y1) || rtIsInf(y1)) {
      y1 = 0.0;
    } else {
      y1 = fmod(y1, 4.294967296E+9);
    }

    tseed = y1 < 0.0 ? (uint32_T)-(int32_T)(uint32_T)-y1 : (uint32_T)y1;
    r = (int32_T)(tseed >> 16U);
    t = (int32_T)(tseed & 32768U);
    tseed = ((((tseed - ((uint32_T)r << 16U)) + t) << 16U) + t) + r;
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else {
      if (tseed > 2147483646U) {
        tseed = 2147483646U;
      }
    }

    y1 = rt_nrand_Upu32_Yd_f_pw_snf(&tseed) * mainV03_56_P.WhiteNoise_StdDev_n +
      mainV03_56_P.WhiteNoise_Mean_d;
    mainV03_56_DW.NextOutput_e[0] = y1;
    mainV03_56_DW.RandSeed_g[0] = tseed;
    y1 = floor(mainV03_56_P.WhiteNoise_Seed_n[1]);
    if (rtIsNaN(y1) || rtIsInf(y1)) {
      y1 = 0.0;
    } else {
      y1 = fmod(y1, 4.294967296E+9);
    }

    tseed = y1 < 0.0 ? (uint32_T)-(int32_T)(uint32_T)-y1 : (uint32_T)y1;
    r = (int32_T)(tseed >> 16U);
    t = (int32_T)(tseed & 32768U);
    tseed = ((((tseed - ((uint32_T)r << 16U)) + t) << 16U) + t) + r;
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else {
      if (tseed > 2147483646U) {
        tseed = 2147483646U;
      }
    }

    y1 = rt_nrand_Upu32_Yd_f_pw_snf(&tseed) * mainV03_56_P.WhiteNoise_StdDev_n +
      mainV03_56_P.WhiteNoise_Mean_d;
    mainV03_56_DW.NextOutput_e[1] = y1;
    mainV03_56_DW.RandSeed_g[1] = tseed;
    y1 = floor(mainV03_56_P.WhiteNoise_Seed_n[2]);
    if (rtIsNaN(y1) || rtIsInf(y1)) {
      y1 = 0.0;
    } else {
      y1 = fmod(y1, 4.294967296E+9);
    }

    tseed = y1 < 0.0 ? (uint32_T)-(int32_T)(uint32_T)-y1 : (uint32_T)y1;
    r = (int32_T)(tseed >> 16U);
    t = (int32_T)(tseed & 32768U);
    tseed = ((((tseed - ((uint32_T)r << 16U)) + t) << 16U) + t) + r;
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else {
      if (tseed > 2147483646U) {
        tseed = 2147483646U;
      }
    }

    y1 = rt_nrand_Upu32_Yd_f_pw_snf(&tseed) * mainV03_56_P.WhiteNoise_StdDev_n +
      mainV03_56_P.WhiteNoise_Mean_d;
    mainV03_56_DW.NextOutput_e[2] = y1;
    mainV03_56_DW.RandSeed_g[2] = tseed;

    /* End of InitializeConditions for RandomNumber: '<S534>/White Noise' */

    /* InitializeConditions for TransferFcn: '<S550>/Transfer Fcn X' */
    mainV03_56_X.TransferFcnX_CSTATE_e[0] = 0.0;

    /* InitializeConditions for TransferFcn: '<S550>/Transfer Fcn Y' */
    mainV03_56_X.TransferFcnY_CSTATE_j[0] = 0.0;

    /* InitializeConditions for TransferFcn: '<S550>/Transfer Fcn Z' */
    mainV03_56_X.TransferFcnZ_CSTATE_l[0] = 0.0;

    /* InitializeConditions for TransferFcn: '<S550>/Transfer Fcn X' */
    mainV03_56_X.TransferFcnX_CSTATE_e[1] = 0.0;

    /* InitializeConditions for TransferFcn: '<S550>/Transfer Fcn Y' */
    mainV03_56_X.TransferFcnY_CSTATE_j[1] = 0.0;

    /* InitializeConditions for TransferFcn: '<S550>/Transfer Fcn Z' */
    mainV03_56_X.TransferFcnZ_CSTATE_l[1] = 0.0;

    /* InitializeConditions for RandomNumber: '<S548>/White Noise' */
    y1 = floor(mainV03_56_P.WhiteNoise_Seed_a[0]);
    if (rtIsNaN(y1) || rtIsInf(y1)) {
      y1 = 0.0;
    } else {
      y1 = fmod(y1, 4.294967296E+9);
    }

    tseed = y1 < 0.0 ? (uint32_T)-(int32_T)(uint32_T)-y1 : (uint32_T)y1;
    r = (int32_T)(tseed >> 16U);
    t = (int32_T)(tseed & 32768U);
    tseed = ((((tseed - ((uint32_T)r << 16U)) + t) << 16U) + t) + r;
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else {
      if (tseed > 2147483646U) {
        tseed = 2147483646U;
      }
    }

    y1 = rt_nrand_Upu32_Yd_f_pw_snf(&tseed) * mainV03_56_P.WhiteNoise_StdDev_l +
      mainV03_56_P.WhiteNoise_Mean_n;
    mainV03_56_DW.NextOutput_b[0] = y1;
    mainV03_56_DW.RandSeed_j[0] = tseed;

    /* InitializeConditions for Integrator: '<S525>/phi theta psi' */
    mainV03_56_X.phithetapsi_CSTATE_e[0] = mainV03_56_P.phithetapsi_IC[0];

    /* InitializeConditions for RandomNumber: '<S548>/White Noise' */
    y1 = floor(mainV03_56_P.WhiteNoise_Seed_a[1]);
    if (rtIsNaN(y1) || rtIsInf(y1)) {
      y1 = 0.0;
    } else {
      y1 = fmod(y1, 4.294967296E+9);
    }

    tseed = y1 < 0.0 ? (uint32_T)-(int32_T)(uint32_T)-y1 : (uint32_T)y1;
    r = (int32_T)(tseed >> 16U);
    t = (int32_T)(tseed & 32768U);
    tseed = ((((tseed - ((uint32_T)r << 16U)) + t) << 16U) + t) + r;
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else {
      if (tseed > 2147483646U) {
        tseed = 2147483646U;
      }
    }

    y1 = rt_nrand_Upu32_Yd_f_pw_snf(&tseed) * mainV03_56_P.WhiteNoise_StdDev_l +
      mainV03_56_P.WhiteNoise_Mean_n;
    mainV03_56_DW.NextOutput_b[1] = y1;
    mainV03_56_DW.RandSeed_j[1] = tseed;

    /* InitializeConditions for Integrator: '<S525>/phi theta psi' */
    mainV03_56_X.phithetapsi_CSTATE_e[1] = mainV03_56_P.phithetapsi_IC[1];

    /* InitializeConditions for RandomNumber: '<S548>/White Noise' */
    y1 = floor(mainV03_56_P.WhiteNoise_Seed_a[2]);
    if (rtIsNaN(y1) || rtIsInf(y1)) {
      y1 = 0.0;
    } else {
      y1 = fmod(y1, 4.294967296E+9);
    }

    tseed = y1 < 0.0 ? (uint32_T)-(int32_T)(uint32_T)-y1 : (uint32_T)y1;
    r = (int32_T)(tseed >> 16U);
    t = (int32_T)(tseed & 32768U);
    tseed = ((((tseed - ((uint32_T)r << 16U)) + t) << 16U) + t) + r;
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else {
      if (tseed > 2147483646U) {
        tseed = 2147483646U;
      }
    }

    y1 = rt_nrand_Upu32_Yd_f_pw_snf(&tseed) * mainV03_56_P.WhiteNoise_StdDev_l +
      mainV03_56_P.WhiteNoise_Mean_n;
    mainV03_56_DW.NextOutput_b[2] = y1;
    mainV03_56_DW.RandSeed_j[2] = tseed;

    /* InitializeConditions for Integrator: '<S525>/phi theta psi' */
    mainV03_56_X.phithetapsi_CSTATE_e[2] = mainV03_56_P.phithetapsi_IC[2];

    /* InitializeConditions for UnitDelay: '<S464>/Unit Delay' */
    mainV03_56_DW.UnitDelay_DSTATE = mainV03_56_P.UnitDelay_InitialCondition_b;

    /* SystemInitialize for Merge: '<S553>/Merge' */
    mainV03_56_B.Merge[0] = mainV03_56_P.Merge_InitialOutput_k;
    mainV03_56_B.Merge[1] = mainV03_56_P.Merge_InitialOutput_k;
    mainV03_56_B.Merge[2] = mainV03_56_P.Merge_InitialOutput_k;
    mainV03_56_B.Merge[3] = mainV03_56_P.Merge_InitialOutput_k;

    /* SystemInitialize for Enabled SubSystem: '<S104>/Convert from geodetic to  spherical coordinates ' */
    /* SystemInitialize for Iterator SubSystem: '<S108>/For Iterator Subsystem' */
    /* InitializeConditions for UnitDelay: '<S155>/Unit Delay1' */
    mainV03_56_DW.UnitDelay1_DSTATE[0] =
      mainV03_56_P.UnitDelay1_InitialCondition_h1;
    mainV03_56_DW.UnitDelay1_DSTATE[1] =
      mainV03_56_P.UnitDelay1_InitialCondition_h1;

    /* End of SystemInitialize for SubSystem: '<S108>/For Iterator Subsystem' */

    /* SystemInitialize for Outport: '<S108>/sp[13]' */
    memcpy(&mainV03_56_B.OutportBufferForsp13[0], &mainV03_56_P.sp13_Y0[0], 13U *
           sizeof(real_T));

    /* SystemInitialize for Outport: '<S108>/cp[13]' */
    memcpy(&mainV03_56_B.OutportBufferForcp13[0], &mainV03_56_P.cp13_Y0[0], 13U *
           sizeof(real_T));

    /* End of SystemInitialize for SubSystem: '<S104>/Convert from geodetic to  spherical coordinates ' */

    /* SystemInitialize for Enabled SubSystem: '<S104>/Convert from geodetic to  spherical coordinates' */
    /* SystemInitialize for Outport: '<S107>/r' */
    mainV03_56_B.sqrt_l = mainV03_56_P.r_Y0;

    /* SystemInitialize for Outport: '<S107>/ct' */
    mainV03_56_B.Product4_l = mainV03_56_P.ct_Y0;

    /* SystemInitialize for Outport: '<S107>/st' */
    mainV03_56_B.sqrt_m = mainV03_56_P.st_Y0;

    /* SystemInitialize for Outport: '<S107>/sa' */
    mainV03_56_B.Product12 = mainV03_56_P.sa_Y0;

    /* SystemInitialize for Outport: '<S107>/ca' */
    mainV03_56_B.Product11_j = mainV03_56_P.ca_Y0;

    /* End of SystemInitialize for SubSystem: '<S104>/Convert from geodetic to  spherical coordinates' */

    /* SystemInitialize for Iterator SubSystem: '<S104>/Compute magnetic vector in spherical coordinates' */
    /* SystemInitialize for Iterator SubSystem: '<S106>/For Iterator Subsystem' */
    /* InitializeConditions for UnitDelay: '<S115>/Unit Delay1' */
    mainV03_56_DW.UnitDelay1_DSTATE_m =
      mainV03_56_P.UnitDelay1_InitialCondition_a;

    /* InitializeConditions for UnitDelay: '<S115>/Unit Delay3' */
    mainV03_56_DW.UnitDelay3_DSTATE = mainV03_56_P.UnitDelay3_InitialCondition;

    /* InitializeConditions for UnitDelay: '<S115>/Unit Delay2' */
    mainV03_56_DW.UnitDelay2_DSTATE_i = mainV03_56_P.UnitDelay2_InitialCondition;

    /* InitializeConditions for UnitDelay: '<S115>/Unit Delay4' */
    mainV03_56_DW.UnitDelay4_DSTATE = mainV03_56_P.UnitDelay4_InitialCondition;

    /* SystemInitialize for Enabled SubSystem: '<S114>/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations' */
    /* SystemInitialize for Merge: '<S116>/Merge1' */
    mainV03_56_B.Merge1_j = mainV03_56_P.Merge1_InitialOutput;

    /* SystemInitialize for Merge: '<S116>/Merge' */
    mainV03_56_B.Merge_im = mainV03_56_P.Merge_InitialOutput;

    /* End of SystemInitialize for SubSystem: '<S114>/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations' */

    /* SystemInitialize for Enabled SubSystem: '<S114>/Time adjust the gauss coefficients' */
    /* InitializeConditions for UnitDelay: '<S117>/Unit Delay' */
    memcpy(&mainV03_56_DW.UnitDelay_DSTATE_f[0],
           &mainV03_56_P.UnitDelay_InitialCondition_e2[0], 169U * sizeof(real_T));

    /* InitializeConditions for UnitDelay: '<S144>/Unit Delay' */
    memcpy(&mainV03_56_DW.UnitDelay_DSTATE_i[0],
           &mainV03_56_P.UnitDelay_InitialCondition_o[0], 169U * sizeof(real_T));

    /* SystemInitialize for Outport: '<S117>/tc[13][13]' */
    memcpy(&mainV03_56_B.Sum2_d[0], &mainV03_56_P.tc1313_Y0[0], 169U * sizeof
           (real_T));

    /* End of SystemInitialize for SubSystem: '<S114>/Time adjust the gauss coefficients' */

    /* SystemInitialize for Enabled SubSystem: '<S114>/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations' */
    /* InitializeConditions for UnitDelay: '<S116>/Unit Delay' */
    memcpy(&mainV03_56_DW.UnitDelay_DSTATE_k[0],
           &mainV03_56_P.UnitDelay_InitialCondition_h[0], 169U * sizeof(real_T));

    /* InitializeConditions for UnitDelay: '<S116>/Unit Delay1' */
    memcpy(&mainV03_56_DW.UnitDelay1_DSTATE_k[0],
           &mainV03_56_P.UnitDelay1_InitialCondition_g[0], 169U * sizeof(real_T));

    /* SystemInitialize for Outport: '<S116>/dp[13][13]' */
    memcpy(&mainV03_56_B.Assignment_i[0], &mainV03_56_P.dp1313_Y0[0], 169U *
           sizeof(real_T));

    /* SystemInitialize for Outport: '<S116>/snorm[169]' */
    memcpy(&mainV03_56_B.Assignment_snorm[0], &mainV03_56_P.snorm169_Y0[0], 169U
           * sizeof(real_T));

    /* End of SystemInitialize for SubSystem: '<S114>/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations' */

    /* SystemInitialize for Enabled SubSystem: '<S115>/Special case - North//South Geographic Pole' */
    /* InitializeConditions for UnitDelay: '<S118>/Unit Delay1' */
    memcpy(&mainV03_56_DW.UnitDelay1_DSTATE_n[0],
           &mainV03_56_P.UnitDelay1_InitialCondition_l[0], 13U * sizeof(real_T));

    /* SystemInitialize for IfAction SubSystem: '<S118>/If Action Subsystem1' */
    /* SystemInitialize for Outport: '<S122>/pp[13]' */
    memcpy(&mainV03_56_B.Assignment2_e[0], &mainV03_56_P.pp13_Y0[0], 13U *
           sizeof(real_T));

    /* End of SystemInitialize for SubSystem: '<S118>/If Action Subsystem1' */

    /* SystemInitialize for IfAction SubSystem: '<S118>/If Action Subsystem2' */
    /* SystemInitialize for Outport: '<S123>/pp[13]' */
    memcpy(&mainV03_56_B.Assignment2_o[0], &mainV03_56_P.pp13_Y0_d[0], 13U *
           sizeof(real_T));

    /* End of SystemInitialize for SubSystem: '<S118>/If Action Subsystem2' */

    /* SystemInitialize for Outport: '<S118>/bpp' */
    mainV03_56_B.Product2_o = mainV03_56_P.bpp_Y0;

    /* End of SystemInitialize for SubSystem: '<S115>/Special case - North//South Geographic Pole' */

    /* SystemInitialize for Outport: '<S114>/bt' */
    mainV03_56_B.Sum1_bu = mainV03_56_P.bt_Y0;

    /* SystemInitialize for Outport: '<S114>/bp' */
    mainV03_56_B.Sum2_p = mainV03_56_P.bp_Y0;

    /* SystemInitialize for Outport: '<S114>/br' */
    mainV03_56_B.Sum3_co = mainV03_56_P.br_Y0;

    /* SystemInitialize for Outport: '<S114>/bpp' */
    mainV03_56_B.Sum5_k = mainV03_56_P.bpp_Y0_f;

    /* End of SystemInitialize for SubSystem: '<S106>/For Iterator Subsystem' */

    /* SystemInitialize for Outport: '<S106>/bt,bp,br,bpp' */
    mainV03_56_B.Sum1_lw[0] = mainV03_56_P.btbpbrbpp_Y0[0];
    mainV03_56_B.Sum1_lw[1] = mainV03_56_P.btbpbrbpp_Y0[1];
    mainV03_56_B.Sum1_lw[2] = mainV03_56_P.btbpbrbpp_Y0[2];
    mainV03_56_B.Sum1_lw[3] = mainV03_56_P.btbpbrbpp_Y0[3];

    /* End of SystemInitialize for SubSystem: '<S104>/Compute magnetic vector in spherical coordinates' */

    /* SystemInitialize for Enabled SubSystem: '<S47>/Distance into gust (x)' */
    /* InitializeConditions for Integrator: '<S50>/Distance into Gust (x) (Limited to gust length d)' */
    mainV03_56_X.DistanceintoGustxLimitedtogustlengthd_CSTATE_a =
      mainV03_56_P.DistanceintoGustxLimitedtogustlengthd_IC;

    /* SystemInitialize for Outport: '<S50>/x' */
    mainV03_56_B.DistanceintoGustxLimitedtogustlengthd = mainV03_56_P.x_Y0;

    /* End of SystemInitialize for SubSystem: '<S47>/Distance into gust (x)' */

    /* SystemInitialize for Enabled SubSystem: '<S47>/Distance into gust (y)' */
    mainV03_56_Distanceintogusty_Init(&mainV03_56_B.Distanceintogusty,
      (P_Distanceintogusty_mainV03_56_T *)&mainV03_56_P.Distanceintogusty,
      &mainV03_56_X.Distanceintogusty);

    /* End of SystemInitialize for SubSystem: '<S47>/Distance into gust (y)' */

    /* SystemInitialize for Enabled SubSystem: '<S47>/Distance into gust (z)' */
    mainV03_56_Distanceintogusty_Init(&mainV03_56_B.Distanceintogustz,
      (P_Distanceintogusty_mainV03_56_T *)&mainV03_56_P.Distanceintogustz,
      &mainV03_56_X.Distanceintogustz);

    /* End of SystemInitialize for SubSystem: '<S47>/Distance into gust (z)' */

    /* SystemInitialize for Enabled SubSystem: '<S55>/Hpgw' */
    /* InitializeConditions for UnitDelay: '<S67>/Unit Delay' */
    mainV03_56_DW.UnitDelay_DSTATE_c[0] =
      mainV03_56_P.UnitDelay_InitialCondition;

    /* SystemInitialize for Outport: '<S67>/pgw' */
    mainV03_56_B.Sum_ax[0] = mainV03_56_P.pgw_Y0;

    /* End of SystemInitialize for SubSystem: '<S55>/Hpgw' */

    /* SystemInitialize for Enabled SubSystem: '<S56>/Hwgw(z)' */
    /* InitializeConditions for UnitDelay: '<S72>/Unit Delay' */
    mainV03_56_DW.UnitDelay_DSTATE_l[0] =
      mainV03_56_P.UnitDelay_InitialCondition_kv;

    /* SystemInitialize for Outport: '<S72>/wgw' */
    mainV03_56_B.Sum_da[0] = mainV03_56_P.wgw_Y0;

    /* End of SystemInitialize for SubSystem: '<S56>/Hwgw(z)' */

    /* SystemInitialize for Enabled SubSystem: '<S55>/Hqgw' */
    /* InitializeConditions for UnitDelay: '<S68>/Unit Delay' */
    mainV03_56_DW.UnitDelay_DSTATE_fj[0] =
      mainV03_56_P.UnitDelay_InitialCondition_f;

    /* InitializeConditions for UnitDelay: '<S68>/Unit Delay1' */
    mainV03_56_DW.UnitDelay1_DSTATE_j[0] =
      mainV03_56_P.UnitDelay1_InitialCondition;

    /* SystemInitialize for Outport: '<S68>/qgw' */
    mainV03_56_B.Sum1_gx[0] = mainV03_56_P.qgw_Y0;

    /* End of SystemInitialize for SubSystem: '<S55>/Hqgw' */

    /* SystemInitialize for Enabled SubSystem: '<S56>/Hvgw(z)' */
    /* InitializeConditions for UnitDelay: '<S71>/Unit Delay' */
    mainV03_56_DW.UnitDelay_DSTATE_o[0] =
      mainV03_56_P.UnitDelay_InitialCondition_e;

    /* SystemInitialize for Outport: '<S71>/vgw' */
    mainV03_56_B.Sum_nt[0] = mainV03_56_P.vgw_Y0;

    /* End of SystemInitialize for SubSystem: '<S56>/Hvgw(z)' */

    /* SystemInitialize for Enabled SubSystem: '<S55>/Hrgw' */
    /* InitializeConditions for UnitDelay: '<S69>/Unit Delay' */
    mainV03_56_DW.UnitDelay_DSTATE_ju[0] =
      mainV03_56_P.UnitDelay_InitialCondition_k;

    /* InitializeConditions for UnitDelay: '<S69>/Unit Delay1' */
    mainV03_56_DW.UnitDelay1_DSTATE_ml[0] =
      mainV03_56_P.UnitDelay1_InitialCondition_h;

    /* SystemInitialize for Outport: '<S69>/rgw' */
    mainV03_56_B.Sum1_p2[0] = mainV03_56_P.rgw_Y0;

    /* End of SystemInitialize for SubSystem: '<S55>/Hrgw' */

    /* SystemInitialize for Enabled SubSystem: '<S56>/Hugw(z)' */
    /* InitializeConditions for UnitDelay: '<S70>/Unit Delay' */
    mainV03_56_DW.UnitDelay_DSTATE_p[0] =
      mainV03_56_P.UnitDelay_InitialCondition_m;

    /* SystemInitialize for Outport: '<S70>/ugw' */
    mainV03_56_B.Sum_eq[0] = mainV03_56_P.ugw_Y0;

    /* End of SystemInitialize for SubSystem: '<S56>/Hugw(z)' */

    /* SystemInitialize for Enabled SubSystem: '<S55>/Hpgw' */
    /* InitializeConditions for UnitDelay: '<S67>/Unit Delay' */
    mainV03_56_DW.UnitDelay_DSTATE_c[1] =
      mainV03_56_P.UnitDelay_InitialCondition;

    /* SystemInitialize for Outport: '<S67>/pgw' */
    mainV03_56_B.Sum_ax[1] = mainV03_56_P.pgw_Y0;

    /* End of SystemInitialize for SubSystem: '<S55>/Hpgw' */

    /* SystemInitialize for Enabled SubSystem: '<S56>/Hwgw(z)' */
    /* InitializeConditions for UnitDelay: '<S72>/Unit Delay' */
    mainV03_56_DW.UnitDelay_DSTATE_l[1] =
      mainV03_56_P.UnitDelay_InitialCondition_kv;

    /* SystemInitialize for Outport: '<S72>/wgw' */
    mainV03_56_B.Sum_da[1] = mainV03_56_P.wgw_Y0;

    /* End of SystemInitialize for SubSystem: '<S56>/Hwgw(z)' */

    /* SystemInitialize for Enabled SubSystem: '<S55>/Hqgw' */
    /* InitializeConditions for UnitDelay: '<S68>/Unit Delay' */
    mainV03_56_DW.UnitDelay_DSTATE_fj[1] =
      mainV03_56_P.UnitDelay_InitialCondition_f;

    /* InitializeConditions for UnitDelay: '<S68>/Unit Delay1' */
    mainV03_56_DW.UnitDelay1_DSTATE_j[1] =
      mainV03_56_P.UnitDelay1_InitialCondition;

    /* SystemInitialize for Outport: '<S68>/qgw' */
    mainV03_56_B.Sum1_gx[1] = mainV03_56_P.qgw_Y0;

    /* End of SystemInitialize for SubSystem: '<S55>/Hqgw' */

    /* SystemInitialize for Enabled SubSystem: '<S56>/Hvgw(z)' */
    /* InitializeConditions for UnitDelay: '<S71>/Unit Delay' */
    mainV03_56_DW.UnitDelay_DSTATE_o[1] =
      mainV03_56_P.UnitDelay_InitialCondition_e;

    /* SystemInitialize for Outport: '<S71>/vgw' */
    mainV03_56_B.Sum_nt[1] = mainV03_56_P.vgw_Y0;

    /* End of SystemInitialize for SubSystem: '<S56>/Hvgw(z)' */

    /* SystemInitialize for Enabled SubSystem: '<S55>/Hrgw' */
    /* InitializeConditions for UnitDelay: '<S69>/Unit Delay' */
    mainV03_56_DW.UnitDelay_DSTATE_ju[1] =
      mainV03_56_P.UnitDelay_InitialCondition_k;

    /* InitializeConditions for UnitDelay: '<S69>/Unit Delay1' */
    mainV03_56_DW.UnitDelay1_DSTATE_ml[1] =
      mainV03_56_P.UnitDelay1_InitialCondition_h;

    /* SystemInitialize for Outport: '<S69>/rgw' */
    mainV03_56_B.Sum1_p2[1] = mainV03_56_P.rgw_Y0;

    /* End of SystemInitialize for SubSystem: '<S55>/Hrgw' */

    /* SystemInitialize for Enabled SubSystem: '<S56>/Hugw(z)' */
    /* InitializeConditions for UnitDelay: '<S70>/Unit Delay' */
    mainV03_56_DW.UnitDelay_DSTATE_p[1] =
      mainV03_56_P.UnitDelay_InitialCondition_m;

    /* SystemInitialize for Outport: '<S70>/ugw' */
    mainV03_56_B.Sum_eq[1] = mainV03_56_P.ugw_Y0;

    /* End of SystemInitialize for SubSystem: '<S56>/Hugw(z)' */

    /* SystemInitialize for Enabled SubSystem: '<S3>/Quadcopter' */
    /* InitializeConditions for Integrator: '<S179>/Filter' */
    mainV03_56_X.Filter_CSTATE_ci = mainV03_56_P.Filter_IC_j;

    /* InitializeConditions for Integrator: '<S179>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_nu = mainV03_56_P.Integrator_IC_i;

    /* InitializeConditions for Integrator: '<S184>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_aw = mainV03_56_P.Integrator_IC_o0;

    /* InitializeConditions for Integrator: '<S184>/Filter' */
    mainV03_56_X.Filter_CSTATE_aw = mainV03_56_P.Filter_IC_eq;

    /* InitializeConditions for Integrator: '<S180>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_ij = mainV03_56_P.Integrator_IC_p;

    /* InitializeConditions for Integrator: '<S180>/Filter' */
    mainV03_56_X.Filter_CSTATE_jg = mainV03_56_P.Filter_IC_b;

    /* InitializeConditions for Integrator: '<S185>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_hz = mainV03_56_P.Integrator_IC_jo;

    /* InitializeConditions for Integrator: '<S185>/Filter' */
    mainV03_56_X.Filter_CSTATE_l3 = mainV03_56_P.Filter_IC_eg;

    /* InitializeConditions for Integrator: '<S183>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_pc = mainV03_56_P.Integrator_IC_k;

    /* InitializeConditions for Integrator: '<S183>/Filter' */
    mainV03_56_X.Filter_CSTATE_c0 = mainV03_56_P.Filter_IC_h;

    /* InitializeConditions for Integrator: '<S191>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_l2 = mainV03_56_P.Integrator_IC_jb;

    /* InitializeConditions for Integrator: '<S191>/Filter' */
    mainV03_56_X.Filter_CSTATE_a5 = mainV03_56_P.Filter_IC_f1;

    /* InitializeConditions for Integrator: '<S187>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_am = mainV03_56_P.Integrator_IC_m;

    /* InitializeConditions for Integrator: '<S187>/Filter' */
    mainV03_56_X.Filter_CSTATE_mg = mainV03_56_P.Filter_IC_cm;

    /* InitializeConditions for Integrator: '<S192>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_dkb = mainV03_56_P.Integrator_IC_jt;

    /* InitializeConditions for Integrator: '<S192>/Filter' */
    mainV03_56_X.Filter_CSTATE_oz = mainV03_56_P.Filter_IC_hb;

    /* InitializeConditions for Integrator: '<S188>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_i5 = mainV03_56_P.Integrator_IC_p2;

    /* InitializeConditions for Integrator: '<S188>/Filter' */
    mainV03_56_X.Filter_CSTATE_px = mainV03_56_P.Filter_IC_or;

    /* InitializeConditions for Integrator: '<S182>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_ng = mainV03_56_P.Integrator_IC_ir;

    /* InitializeConditions for Integrator: '<S182>/Filter' */
    mainV03_56_X.Filter_CSTATE_d = mainV03_56_P.Filter_IC_p;

    /* InitializeConditions for Integrator: '<S190>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_cp = mainV03_56_P.Integrator_IC_fd;

    /* InitializeConditions for Integrator: '<S190>/Filter' */
    mainV03_56_X.Filter_CSTATE_gi = mainV03_56_P.Filter_IC_ao;

    /* InitializeConditions for Integrator: '<S186>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_e1 = mainV03_56_P.Integrator_IC_ev;

    /* InitializeConditions for Integrator: '<S186>/Filter' */
    mainV03_56_X.Filter_CSTATE_pc = mainV03_56_P.Filter_IC_i;

    /* InitializeConditions for Integrator: '<S189>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_gw = mainV03_56_P.Integrator_IC_ph;

    /* InitializeConditions for Integrator: '<S189>/Filter' */
    mainV03_56_X.Filter_CSTATE_py = mainV03_56_P.Filter_IC_pe;

    /* InitializeConditions for Integrator: '<S181>/Filter' */
    mainV03_56_X.Filter_CSTATE_ib = mainV03_56_P.Filter_IC_g;

    /* InitializeConditions for Integrator: '<S181>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_l4 = mainV03_56_P.Integrator_IC_nq;

    /* SystemInitialize for Outport: '<S165>/actuators' */
    mainV03_56_B.actuators_h = mainV03_56_P.actuators_Y0_n;

    /* SystemInitialize for Outport: '<S165>/Throttle' */
    mainV03_56_B.ManualSwitch3_i = mainV03_56_P.Throttle_Y0_b;

    /* End of SystemInitialize for SubSystem: '<S3>/Quadcopter' */

    /* SystemInitialize for Enabled SubSystem: '<S3>/Quadcopter --> Fixed-Wing' */
    /* InitializeConditions for Integrator: '<S201>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_n1 = mainV03_56_P.Integrator_IC_kg;

    /* InitializeConditions for Integrator: '<S201>/Filter' */
    mainV03_56_X.Filter_CSTATE_nt = mainV03_56_P.Filter_IC_l;

    /* InitializeConditions for Integrator: '<S207>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_e = mainV03_56_P.Integrator_IC_ig;

    /* InitializeConditions for Integrator: '<S207>/Filter' */
    mainV03_56_X.Filter_CSTATE_m = mainV03_56_P.Filter_IC_pt;

    /* InitializeConditions for Integrator: '<S206>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_il = mainV03_56_P.Integrator_IC_nu;

    /* InitializeConditions for Integrator: '<S206>/Filter' */
    mainV03_56_X.Filter_CSTATE_i = mainV03_56_P.Filter_IC_ms;

    /* InitializeConditions for Integrator: '<S209>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_mj = mainV03_56_P.Integrator_IC_l;

    /* InitializeConditions for Integrator: '<S209>/Filter' */
    mainV03_56_X.Filter_CSTATE_kd = mainV03_56_P.Filter_IC_em;

    /* InitializeConditions for Integrator: '<S208>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_dk = mainV03_56_P.Integrator_IC_np;

    /* InitializeConditions for Integrator: '<S208>/Filter' */
    mainV03_56_X.Filter_CSTATE_ii = mainV03_56_P.Filter_IC_gj;

    /* InitializeConditions for Integrator: '<S211>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_l = mainV03_56_P.Integrator_IC_ot;

    /* InitializeConditions for Integrator: '<S211>/Filter' */
    mainV03_56_X.Filter_CSTATE_c = mainV03_56_P.Filter_IC_lx;

    /* InitializeConditions for Integrator: '<S210>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_p0 = mainV03_56_P.Integrator_IC_b;

    /* InitializeConditions for Integrator: '<S210>/Filter' */
    mainV03_56_X.Filter_CSTATE_a = mainV03_56_P.Filter_IC_gm;

    /* InitializeConditions for Integrator: '<S200>/Filter' */
    mainV03_56_X.Filter_CSTATE_mr = mainV03_56_P.Filter_IC_is;

    /* InitializeConditions for Integrator: '<S200>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_mx = mainV03_56_P.Integrator_IC_br;

    /* InitializeConditions for Integrator: '<S212>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_g = mainV03_56_P.Integrator_IC_n5;

    /* InitializeConditions for Integrator: '<S212>/Filter' */
    mainV03_56_X.Filter_CSTATE_o = mainV03_56_P.Filter_IC_gn;

    /* InitializeConditions for Integrator: '<S202>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_o = mainV03_56_P.Integrator_IC_fx;

    /* InitializeConditions for Integrator: '<S202>/Filter' */
    mainV03_56_X.Filter_CSTATE_mk = mainV03_56_P.Filter_IC_d;

    /* InitializeConditions for Integrator: '<S213>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_io = mainV03_56_P.Integrator_IC_a;

    /* InitializeConditions for Integrator: '<S213>/Filter' */
    mainV03_56_X.Filter_CSTATE_p = mainV03_56_P.Filter_IC_mu;

    /* InitializeConditions for Integrator: '<S205>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_m3 = mainV03_56_P.Integrator_IC_ol;

    /* InitializeConditions for Integrator: '<S205>/Filter' */
    mainV03_56_X.Filter_CSTATE_kdq = mainV03_56_P.Filter_IC_jf;

    /* InitializeConditions for Integrator: '<S219>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_gq = mainV03_56_P.Integrator_IC_c;

    /* InitializeConditions for Integrator: '<S219>/Filter' */
    mainV03_56_X.Filter_CSTATE_az = mainV03_56_P.Filter_IC_og;

    /* InitializeConditions for Integrator: '<S215>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_nb = mainV03_56_P.Integrator_IC_l5;

    /* InitializeConditions for Integrator: '<S215>/Filter' */
    mainV03_56_X.Filter_CSTATE_f = mainV03_56_P.Filter_IC_gmg;

    /* InitializeConditions for Integrator: '<S220>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_j = mainV03_56_P.Integrator_IC_mg;

    /* InitializeConditions for Integrator: '<S220>/Filter' */
    mainV03_56_X.Filter_CSTATE_ok = mainV03_56_P.Filter_IC_n;

    /* InitializeConditions for Integrator: '<S216>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_h = mainV03_56_P.Integrator_IC_b2;

    /* InitializeConditions for Integrator: '<S216>/Filter' */
    mainV03_56_X.Filter_CSTATE_or = mainV03_56_P.Filter_IC_fw;

    /* InitializeConditions for Integrator: '<S204>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_b = mainV03_56_P.Integrator_IC_e1;

    /* InitializeConditions for Integrator: '<S204>/Filter' */
    mainV03_56_X.Filter_CSTATE_j = mainV03_56_P.Filter_IC_cy;

    /* InitializeConditions for Integrator: '<S218>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_e5 = mainV03_56_P.Integrator_IC_ob;

    /* InitializeConditions for Integrator: '<S218>/Filter' */
    mainV03_56_X.Filter_CSTATE_mkr = mainV03_56_P.Filter_IC_hf;

    /* InitializeConditions for Integrator: '<S214>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_ni = mainV03_56_P.Integrator_IC_po;

    /* InitializeConditions for Integrator: '<S214>/Filter' */
    mainV03_56_X.Filter_CSTATE_l = mainV03_56_P.Filter_IC_oz;

    /* InitializeConditions for Integrator: '<S217>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_c = mainV03_56_P.Integrator_IC_ey;

    /* InitializeConditions for Integrator: '<S217>/Filter' */
    mainV03_56_X.Filter_CSTATE_j2 = mainV03_56_P.Filter_IC_od;

    /* InitializeConditions for Integrator: '<S203>/Filter' */
    mainV03_56_X.Filter_CSTATE_br = mainV03_56_P.Filter_IC_cm1;

    /* InitializeConditions for Integrator: '<S203>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_a = mainV03_56_P.Integrator_IC_pf;

    /* SystemInitialize for Outport: '<S166>/actuators' */
    mainV03_56_B.BusCreator2 = mainV03_56_P.actuators_Y0_h;

    /* SystemInitialize for Outport: '<S166>/Throttle' */
    mainV03_56_B.ManualSwitch3 = mainV03_56_P.Throttle_Y0_e;

    /* End of SystemInitialize for SubSystem: '<S3>/Quadcopter --> Fixed-Wing' */

    /* SystemInitialize for Enabled SubSystem: '<S3>/Fixed-Wing Climb' */
    /* InitializeConditions for Integrator: '<S176>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_dg = mainV03_56_P.Integrator_IC;

    /* InitializeConditions for Integrator: '<S176>/Filter' */
    mainV03_56_X.Filter_CSTATE_k1 = mainV03_56_P.Filter_IC;

    /* InitializeConditions for Integrator: '<S169>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_ke = mainV03_56_P.Integrator_IC_o;

    /* InitializeConditions for Integrator: '<S169>/Filter' */
    mainV03_56_X.Filter_CSTATE_ba = mainV03_56_P.Filter_IC_m;

    /* InitializeConditions for Integrator: '<S171>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_ih = mainV03_56_P.Integrator_IC_d;

    /* InitializeConditions for Integrator: '<S171>/Filter' */
    mainV03_56_X.Filter_CSTATE_lp = mainV03_56_P.Filter_IC_a;

    /* InitializeConditions for Integrator: '<S170>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_nm = mainV03_56_P.Integrator_IC_j;

    /* InitializeConditions for Integrator: '<S170>/Filter' */
    mainV03_56_X.Filter_CSTATE_dt = mainV03_56_P.Filter_IC_o;

    /* InitializeConditions for Integrator: '<S173>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_og = mainV03_56_P.Integrator_IC_e;

    /* InitializeConditions for Integrator: '<S173>/Filter' */
    mainV03_56_X.Filter_CSTATE_iv = mainV03_56_P.Filter_IC_ar;

    /* InitializeConditions for Integrator: '<S172>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_ht = mainV03_56_P.Integrator_IC_n;

    /* InitializeConditions for Integrator: '<S172>/Filter' */
    mainV03_56_X.Filter_CSTATE_k5 = mainV03_56_P.Filter_IC_f;

    /* InitializeConditions for Integrator: '<S175>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_cw = mainV03_56_P.Integrator_IC_f;

    /* InitializeConditions for Integrator: '<S175>/Filter' */
    mainV03_56_X.Filter_CSTATE_lo = mainV03_56_P.Filter_IC_e;

    /* InitializeConditions for Integrator: '<S174>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_ki = mainV03_56_P.Integrator_IC_fq;

    /* InitializeConditions for Integrator: '<S174>/Filter' */
    mainV03_56_X.Filter_CSTATE_kc = mainV03_56_P.Filter_IC_c;

    /* SystemInitialize for Outport: '<S164>/actuators' */
    mainV03_56_B.ManualSwitch1_n = mainV03_56_P.actuators_Y0;

    /* SystemInitialize for Outport: '<S164>/Throttle' */
    mainV03_56_B.ManualSwitch5_j = mainV03_56_P.Throttle_Y0;

    /* End of SystemInitialize for SubSystem: '<S3>/Fixed-Wing Climb' */

    /* SystemInitialize for Enabled SubSystem: '<S3>/Fixed-Wing - Cruise' */
    /* SystemInitialize for BusCreator: '<S163>/Bus Creator1' */
    mainV03_56_B.Throttle_m.Throttle1 = mainV03_56_B.Throttle1;
    mainV03_56_B.Throttle_m.Throttle2 = mainV03_56_B.Throttle2_g;
    mainV03_56_B.Throttle_m.Throttle3 = mainV03_56_B.Throttle3_jz;
    mainV03_56_B.Throttle_m.Throttle4 = mainV03_56_B.Throttle4_j;
    mainV03_56_B.Throttle_m.Throttle5 = mainV03_56_B.Throttle5_a;

    /* SystemInitialize for Outport: '<S163>/Actuators' */
    mainV03_56_B.Actuators = mainV03_56_P.Actuators_Y0;

    /* End of SystemInitialize for SubSystem: '<S3>/Fixed-Wing - Cruise' */

    /* SystemInitialize for Enabled SubSystem: '<S3>/TakeOff' */
    /* InitializeConditions for Integrator: '<S233>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_p = mainV03_56_P.Integrator_IC_du;

    /* InitializeConditions for Integrator: '<S233>/Filter' */
    mainV03_56_X.Filter_CSTATE = mainV03_56_P.Filter_IC_eb;

    /* InitializeConditions for Integrator: '<S235>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_i = mainV03_56_P.Integrator_IC_p2w;

    /* InitializeConditions for Integrator: '<S235>/Filter' */
    mainV03_56_X.Filter_CSTATE_e = mainV03_56_P.Filter_IC_lf;

    /* InitializeConditions for Integrator: '<S234>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_n = mainV03_56_P.Integrator_IC_js;

    /* InitializeConditions for Integrator: '<S234>/Filter' */
    mainV03_56_X.Filter_CSTATE_b = mainV03_56_P.Filter_IC_iz;

    /* InitializeConditions for Integrator: '<S237>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_d = mainV03_56_P.Integrator_IC_bm;

    /* InitializeConditions for Integrator: '<S237>/Filter' */
    mainV03_56_X.Filter_CSTATE_g = mainV03_56_P.Filter_IC_ib;

    /* InitializeConditions for Integrator: '<S236>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_k = mainV03_56_P.Integrator_IC_ds;

    /* InitializeConditions for Integrator: '<S236>/Filter' */
    mainV03_56_X.Filter_CSTATE_gp = mainV03_56_P.Filter_IC_ck;

    /* InitializeConditions for Integrator: '<S239>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_m = mainV03_56_P.Integrator_IC_h;

    /* InitializeConditions for Integrator: '<S239>/Filter' */
    mainV03_56_X.Filter_CSTATE_k = mainV03_56_P.Filter_IC_oo;

    /* InitializeConditions for Integrator: '<S238>/Integrator' */
    mainV03_56_X.Integrator_CSTATE_py = mainV03_56_P.Integrator_IC_ep;

    /* InitializeConditions for Integrator: '<S238>/Filter' */
    mainV03_56_X.Filter_CSTATE_n = mainV03_56_P.Filter_IC_pk;

    /* SystemInitialize for Outport: '<S168>/actuators' */
    mainV03_56_B.ManualSwitch1 = mainV03_56_P.actuators_Y0_m;

    /* SystemInitialize for Outport: '<S168>/Throttle' */
    mainV03_56_B.ManualSwitch5 = mainV03_56_P.Throttle_Y0_m;

    /* End of SystemInitialize for SubSystem: '<S3>/TakeOff' */
  }
}

/* Model terminate function */
void mainV03_56_terminate(void)
{
  /* Level2 S-Function Block: '<S5>/Joystick Input' (joyinput) */
  {
    SimStruct *rts = mainV03_56_M->childSfunctions[0];
    sfcnTerminate(rts);
  }
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/

/* Solver interface called by GRT_Main */
#ifndef USE_GENERATED_SOLVER

void rt_ODECreateIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEDestroyIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEUpdateContinuousStates(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

#endif

void MdlOutputs(int_T tid)
{
  mainV03_56_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  mainV03_56_update();
  UNUSED_PARAMETER(tid);
}

void MdlInitializeSizes(void)
{
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
}

void MdlStart(void)
{
  mainV03_56_initialize();
}

void MdlTerminate(void)
{
  mainV03_56_terminate();
}

/* Registration function */
RT_MODEL_mainV03_56_T *mainV03_56(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  mainV03_56_P.ProportionalAltitude_LowerSaturationLimit = rtMinusInf;
  mainV03_56_P.ProportionalAltitude_LowerSaturationLimit_c = rtMinusInf;
  mainV03_56_P.uftinf_UpperSat = rtInf;
  mainV03_56_P.Saturation_UpperSat_a[0] = rtInf;
  mainV03_56_P.Saturation_UpperSat_a[1] = rtInf;
  mainV03_56_P.Saturation_UpperSat_a[2] = rtInf;
  mainV03_56_P.Saturation_LowerSat_o[0] = rtMinusInf;
  mainV03_56_P.Saturation_LowerSat_o[1] = rtMinusInf;
  mainV03_56_P.Saturation_LowerSat_o[2] = rtMinusInf;
  mainV03_56_P.Saturation_UpperSat_c[0] = rtInf;
  mainV03_56_P.Saturation_UpperSat_c[1] = rtInf;
  mainV03_56_P.Saturation_UpperSat_c[2] = rtInf;
  mainV03_56_P.Saturation_LowerSat_m[0] = rtMinusInf;
  mainV03_56_P.Saturation_LowerSat_m[1] = rtMinusInf;
  mainV03_56_P.Saturation_LowerSat_m[2] = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)mainV03_56_M, 0,
                sizeof(RT_MODEL_mainV03_56_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&mainV03_56_M->solverInfo,
                          &mainV03_56_M->Timing.simTimeStep);
    rtsiSetTPtr(&mainV03_56_M->solverInfo, &rtmGetTPtr(mainV03_56_M));
    rtsiSetStepSizePtr(&mainV03_56_M->solverInfo,
                       &mainV03_56_M->Timing.stepSize0);
    rtsiSetdXPtr(&mainV03_56_M->solverInfo, &mainV03_56_M->derivs);
    rtsiSetContStatesPtr(&mainV03_56_M->solverInfo, (real_T **)
                         &mainV03_56_M->contStates);
    rtsiSetNumContStatesPtr(&mainV03_56_M->solverInfo,
      &mainV03_56_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&mainV03_56_M->solverInfo,
      &mainV03_56_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&mainV03_56_M->solverInfo,
      &mainV03_56_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&mainV03_56_M->solverInfo,
      &mainV03_56_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&mainV03_56_M->solverInfo, (&rtmGetErrorStatus
      (mainV03_56_M)));
    rtsiSetRTModelPtr(&mainV03_56_M->solverInfo, mainV03_56_M);
  }

  rtsiSetSimTimeStep(&mainV03_56_M->solverInfo, MAJOR_TIME_STEP);
  mainV03_56_M->intgData.x0 = mainV03_56_M->odeX0;
  mainV03_56_M->intgData.f0 = mainV03_56_M->odeF0;
  mainV03_56_M->intgData.x1start = mainV03_56_M->odeX1START;
  mainV03_56_M->intgData.f1 = mainV03_56_M->odeF1;
  mainV03_56_M->intgData.Delta = mainV03_56_M->odeDELTA;
  mainV03_56_M->intgData.E = mainV03_56_M->odeE;
  mainV03_56_M->intgData.fac = mainV03_56_M->odeFAC;

  /* initialize */
  {
    int_T i;
    real_T *f = mainV03_56_M->intgData.fac;
    for (i = 0; i < (int_T)(sizeof(mainV03_56_M->odeFAC)/sizeof(real_T)); i++) {
      f[i] = 1.5e-8;
    }
  }

  mainV03_56_M->intgData.DFDX = mainV03_56_M->odeDFDX;
  mainV03_56_M->intgData.W = mainV03_56_M->odeW;
  mainV03_56_M->intgData.pivots = mainV03_56_M->odePIVOTS;
  mainV03_56_M->intgData.xtmp = mainV03_56_M->odeXTMP;
  mainV03_56_M->intgData.ztmp = mainV03_56_M->odeZTMP;
  mainV03_56_M->intgData.isFirstStep = true;
  rtsiSetSolverExtrapolationOrder(&mainV03_56_M->solverInfo, 4);
  rtsiSetSolverNumberNewtonIterations(&mainV03_56_M->solverInfo, 1);
  mainV03_56_M->contStates = ((real_T *) &mainV03_56_X);
  rtsiSetSolverData(&mainV03_56_M->solverInfo, (void *)&mainV03_56_M->intgData);
  rtsiSetSolverName(&mainV03_56_M->solverInfo,"ode14x");
  mainV03_56_M->solverInfoPtr = (&mainV03_56_M->solverInfo);

  /* Initialize timing info */
  {
    int_T *mdlTsMap = mainV03_56_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    mdlTsMap[2] = 2;
    mainV03_56_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    mainV03_56_M->Timing.sampleTimes = (&mainV03_56_M->Timing.sampleTimesArray[0]);
    mainV03_56_M->Timing.offsetTimes = (&mainV03_56_M->Timing.offsetTimesArray[0]);

    /* task periods */
    mainV03_56_M->Timing.sampleTimes[0] = (0.0);
    mainV03_56_M->Timing.sampleTimes[1] = (0.02);
    mainV03_56_M->Timing.sampleTimes[2] = (0.1);

    /* task offsets */
    mainV03_56_M->Timing.offsetTimes[0] = (0.0);
    mainV03_56_M->Timing.offsetTimes[1] = (0.0);
    mainV03_56_M->Timing.offsetTimes[2] = (0.0);
  }

  rtmSetTPtr(mainV03_56_M, &mainV03_56_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = mainV03_56_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    mdlSampleHits[2] = 1;
    mainV03_56_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(mainV03_56_M, -1);
  mainV03_56_M->Timing.stepSize0 = 0.02;
  mainV03_56_M->Timing.stepSize1 = 0.02;
  mainV03_56_M->Timing.stepSize2 = 0.1;

  /* Setup for data logging */
  {
    static RTWLogInfo rt_DataLoggingInfo;
    rt_DataLoggingInfo.loggingInterval = NULL;
    mainV03_56_M->rtwLogInfo = &rt_DataLoggingInfo;
  }

  /* Setup for data logging */
  {
    rtliSetLogXSignalInfo(mainV03_56_M->rtwLogInfo, (NULL));
    rtliSetLogXSignalPtrs(mainV03_56_M->rtwLogInfo, (NULL));
    rtliSetLogT(mainV03_56_M->rtwLogInfo, "tout");
    rtliSetLogX(mainV03_56_M->rtwLogInfo, "");
    rtliSetLogXFinal(mainV03_56_M->rtwLogInfo, "");
    rtliSetLogVarNameModifier(mainV03_56_M->rtwLogInfo, "rt_");
    rtliSetLogFormat(mainV03_56_M->rtwLogInfo, 4);
    rtliSetLogMaxRows(mainV03_56_M->rtwLogInfo, 0);
    rtliSetLogDecimation(mainV03_56_M->rtwLogInfo, 1);
    rtliSetLogY(mainV03_56_M->rtwLogInfo, "");
    rtliSetLogYSignalInfo(mainV03_56_M->rtwLogInfo, (NULL));
    rtliSetLogYSignalPtrs(mainV03_56_M->rtwLogInfo, (NULL));
  }

  mainV03_56_M->solverInfoPtr = (&mainV03_56_M->solverInfo);
  mainV03_56_M->Timing.stepSize = (0.02);
  rtsiSetFixedStepSize(&mainV03_56_M->solverInfo, 0.02);
  rtsiSetSolverMode(&mainV03_56_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  mainV03_56_M->blockIO = ((void *) &mainV03_56_B);
  (void) memset(((void *) &mainV03_56_B), 0,
                sizeof(B_mainV03_56_T));

  /* parameters */
  mainV03_56_M->defaultParam = ((real_T *)&mainV03_56_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &mainV03_56_X;
    mainV03_56_M->contStates = (x);
    (void) memset((void *)&mainV03_56_X, 0,
                  sizeof(X_mainV03_56_T));
  }

  /* states (dwork) */
  mainV03_56_M->dwork = ((void *) &mainV03_56_DW);
  (void) memset((void *)&mainV03_56_DW, 0,
                sizeof(DW_mainV03_56_T));

  /* child S-Function registration */
  {
    RTWSfcnInfo *sfcnInfo = &mainV03_56_M->NonInlinedSFcns.sfcnInfo;
    mainV03_56_M->sfcnInfo = (sfcnInfo);
    rtssSetErrorStatusPtr(sfcnInfo, (&rtmGetErrorStatus(mainV03_56_M)));
    rtssSetNumRootSampTimesPtr(sfcnInfo, &mainV03_56_M->Sizes.numSampTimes);
    mainV03_56_M->NonInlinedSFcns.taskTimePtrs[0] = &(rtmGetTPtr(mainV03_56_M)[0]);
    mainV03_56_M->NonInlinedSFcns.taskTimePtrs[1] = &(rtmGetTPtr(mainV03_56_M)[1]);
    mainV03_56_M->NonInlinedSFcns.taskTimePtrs[2] = &(rtmGetTPtr(mainV03_56_M)[2]);
    rtssSetTPtrPtr(sfcnInfo,mainV03_56_M->NonInlinedSFcns.taskTimePtrs);
    rtssSetTStartPtr(sfcnInfo, &rtmGetTStart(mainV03_56_M));
    rtssSetTFinalPtr(sfcnInfo, &rtmGetTFinal(mainV03_56_M));
    rtssSetTimeOfLastOutputPtr(sfcnInfo, &rtmGetTimeOfLastOutput(mainV03_56_M));
    rtssSetStepSizePtr(sfcnInfo, &mainV03_56_M->Timing.stepSize);
    rtssSetStopRequestedPtr(sfcnInfo, &rtmGetStopRequested(mainV03_56_M));
    rtssSetDerivCacheNeedsResetPtr(sfcnInfo, &mainV03_56_M->derivCacheNeedsReset);
    rtssSetZCCacheNeedsResetPtr(sfcnInfo, &mainV03_56_M->zCCacheNeedsReset);
    rtssSetBlkStateChangePtr(sfcnInfo, &mainV03_56_M->blkStateChange);
    rtssSetSampleHitsPtr(sfcnInfo, &mainV03_56_M->Timing.sampleHits);
    rtssSetPerTaskSampleHitsPtr(sfcnInfo,
      &mainV03_56_M->Timing.perTaskSampleHits);
    rtssSetSimModePtr(sfcnInfo, &mainV03_56_M->simMode);
    rtssSetSolverInfoPtr(sfcnInfo, &mainV03_56_M->solverInfoPtr);
  }

  mainV03_56_M->Sizes.numSFcns = (1);

  /* register each child */
  {
    (void) memset((void *)&mainV03_56_M->NonInlinedSFcns.childSFunctions[0], 0,
                  1*sizeof(SimStruct));
    mainV03_56_M->childSfunctions =
      (&mainV03_56_M->NonInlinedSFcns.childSFunctionPtrs[0]);
    mainV03_56_M->childSfunctions[0] =
      (&mainV03_56_M->NonInlinedSFcns.childSFunctions[0]);

    /* Level2 S-Function Block: mainV03_56/<S5>/Joystick Input (joyinput) */
    {
      SimStruct *rts = mainV03_56_M->childSfunctions[0];

      /* timing info */
      time_T *sfcnPeriod = mainV03_56_M->NonInlinedSFcns.Sfcn0.sfcnPeriod;
      time_T *sfcnOffset = mainV03_56_M->NonInlinedSFcns.Sfcn0.sfcnOffset;
      int_T *sfcnTsMap = mainV03_56_M->NonInlinedSFcns.Sfcn0.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts, &mainV03_56_M->NonInlinedSFcns.blkInfo2[0]);
      }

      ssSetRTWSfcnInfo(rts, mainV03_56_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts, &mainV03_56_M->NonInlinedSFcns.methods2[0]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts, &mainV03_56_M->NonInlinedSFcns.methods3[0]);
      }

      /* Allocate memory of model methods 4 */
      {
        ssSetModelMethods4(rts, &mainV03_56_M->NonInlinedSFcns.methods4[0]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts, &mainV03_56_M->NonInlinedSFcns.statesInfo2[0]);
        ssSetPeriodicStatesInfo(rts,
          &mainV03_56_M->NonInlinedSFcns.periodicStatesInfo[0]);
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &mainV03_56_M->NonInlinedSFcns.Sfcn0.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 3);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 4);
          ssSetOutputPortSignal(rts, 0, ((real_T *) mainV03_56_B.axes));
        }

        /* port 1 */
        {
          _ssSetOutputPortNumDimensions(rts, 1, 1);
          ssSetOutputPortWidth(rts, 1, 12);
          ssSetOutputPortSignal(rts, 1, ((boolean_T *) mainV03_56_B.buttons));
        }

        /* port 2 */
        {
          _ssSetOutputPortNumDimensions(rts, 2, 1);
          ssSetOutputPortWidth(rts, 2, 1);
          ssSetOutputPortSignal(rts, 2, ((real_T *) &mainV03_56_B.POV));
        }
      }

      /* path info */
      ssSetModelName(rts, "Joystick Input");
      ssSetPath(rts, "mainV03_56/Pilot Commands/Joystick Input");
      ssSetRTModel(rts,mainV03_56_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &mainV03_56_M->NonInlinedSFcns.Sfcn0.params;
        ssSetSFcnParamsCount(rts, 3);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)mainV03_56_P.JoystickInput_P1_Size);
        ssSetSFcnParam(rts, 1, (mxArray*)mainV03_56_P.JoystickInput_P2_Size);
        ssSetSFcnParam(rts, 2, (mxArray*)mainV03_56_P.JoystickInput_P3_Size);
      }

      /* work vectors */
      ssSetIWork(rts, (int_T *) &mainV03_56_DW.JoystickInput_IWORK[0]);
      ssSetPWork(rts, (void **) &mainV03_56_DW.JoystickInput_PWORK[0]);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &mainV03_56_M->NonInlinedSFcns.Sfcn0.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &mainV03_56_M->NonInlinedSFcns.Sfcn0.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 2);

        /* IWORK */
        ssSetDWorkWidth(rts, 0, 5);
        ssSetDWorkDataType(rts, 0,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &mainV03_56_DW.JoystickInput_IWORK[0]);

        /* PWORK */
        ssSetDWorkWidth(rts, 1, 3);
        ssSetDWorkDataType(rts, 1,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &mainV03_56_DW.JoystickInput_PWORK[0]);
      }

      /* registration */
      joyinput(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 0.02);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortConnected(rts, 1, 1);
      _ssSetOutputPortConnected(rts, 2, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);
      _ssSetOutputPortBeingMerged(rts, 1, 0);
      _ssSetOutputPortBeingMerged(rts, 2, 0);

      /* Update the BufferDstPort flags for each input port */
    }
  }

  /* Initialize Sizes */
  mainV03_56_M->Sizes.numContStates = (140);/* Number of continuous states */
  mainV03_56_M->Sizes.numPeriodicContStates = (0);/* Number of periodic continuous states */
  mainV03_56_M->Sizes.numY = (0);      /* Number of model outputs */
  mainV03_56_M->Sizes.numU = (0);      /* Number of model inputs */
  mainV03_56_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  mainV03_56_M->Sizes.numSampTimes = (3);/* Number of sample times */
  mainV03_56_M->Sizes.numBlocks = (2498);/* Number of blocks */
  mainV03_56_M->Sizes.numBlockIO = (1477);/* Number of block outputs */
  mainV03_56_M->Sizes.numBlockPrms = (2931850);/* Sum of parameter "widths" */
  return mainV03_56_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/

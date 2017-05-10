/*
 * mainV03_56.h
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

#ifndef RTW_HEADER_mainV03_56_h_
#define RTW_HEADER_mainV03_56_h_
#include <math.h>
#include <stddef.h>
#include <string.h>
#include <float.h>
#ifndef mainV03_56_COMMON_INCLUDES_
# define mainV03_56_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "zero_crossing_types.h"
#include "simstruc.h"
#include "fixedpoint.h"
#include "rt_logging.h"
#endif                                 /* mainV03_56_COMMON_INCLUDES_ */

#include "mainV03_56_types.h"

/* Shared type includes */
#include "multiword_types.h"
#include "rtGetInf.h"
#include "rt_nonfinite.h"
#include "rt_assert.h"
#include "rt_defines.h"
#include "math.h"
#include "rt_matrixlib.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetBlkStateChangeFlag
# define rtmGetBlkStateChangeFlag(rtm) ((rtm)->blkStateChange)
#endif

#ifndef rtmSetBlkStateChangeFlag
# define rtmSetBlkStateChangeFlag(rtm, val) ((rtm)->blkStateChange = (val))
#endif

#ifndef rtmGetBlockIO
# define rtmGetBlockIO(rtm)            ((rtm)->blockIO)
#endif

#ifndef rtmSetBlockIO
# define rtmSetBlockIO(rtm, val)       ((rtm)->blockIO = (val))
#endif

#ifndef rtmGetChecksums
# define rtmGetChecksums(rtm)          ((rtm)->Sizes.checksums)
#endif

#ifndef rtmSetChecksums
# define rtmSetChecksums(rtm, val)     ((rtm)->Sizes.checksums = (val))
#endif

#ifndef rtmGetConstBlockIO
# define rtmGetConstBlockIO(rtm)       ((rtm)->constBlockIO)
#endif

#ifndef rtmSetConstBlockIO
# define rtmSetConstBlockIO(rtm, val)  ((rtm)->constBlockIO = (val))
#endif

#ifndef rtmGetContStateDisabled
# define rtmGetContStateDisabled(rtm)  ((rtm)->contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
# define rtmSetContStateDisabled(rtm, val) ((rtm)->contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
# define rtmGetContStates(rtm)         ((rtm)->contStates)
#endif

#ifndef rtmSetContStates
# define rtmSetContStates(rtm, val)    ((rtm)->contStates = (val))
#endif

#ifndef rtmGetDataMapInfo
# define rtmGetDataMapInfo(rtm)        ()
#endif

#ifndef rtmSetDataMapInfo
# define rtmSetDataMapInfo(rtm, val)   ()
#endif

#ifndef rtmGetDefaultParam
# define rtmGetDefaultParam(rtm)       ((rtm)->defaultParam)
#endif

#ifndef rtmSetDefaultParam
# define rtmSetDefaultParam(rtm, val)  ((rtm)->defaultParam = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
# define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
# define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetDirectFeedThrough
# define rtmGetDirectFeedThrough(rtm)  ((rtm)->Sizes.sysDirFeedThru)
#endif

#ifndef rtmSetDirectFeedThrough
# define rtmSetDirectFeedThrough(rtm, val) ((rtm)->Sizes.sysDirFeedThru = (val))
#endif

#ifndef rtmGetErrorStatusFlag
# define rtmGetErrorStatusFlag(rtm)    ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatusFlag
# define rtmSetErrorStatusFlag(rtm, val) ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetFinalTime
# define rtmGetFinalTime(rtm)          ((rtm)->Timing.tFinal)
#endif

#ifndef rtmSetFinalTime
# define rtmSetFinalTime(rtm, val)     ((rtm)->Timing.tFinal = (val))
#endif

#ifndef rtmGetFirstInitCondFlag
# define rtmGetFirstInitCondFlag(rtm)  ()
#endif

#ifndef rtmSetFirstInitCondFlag
# define rtmSetFirstInitCondFlag(rtm, val) ()
#endif

#ifndef rtmGetIntgData
# define rtmGetIntgData(rtm)           ((rtm)->intgData)
#endif

#ifndef rtmSetIntgData
# define rtmSetIntgData(rtm, val)      ((rtm)->intgData = (val))
#endif

#ifndef rtmGetMdlRefGlobalTID
# define rtmGetMdlRefGlobalTID(rtm)    ()
#endif

#ifndef rtmSetMdlRefGlobalTID
# define rtmSetMdlRefGlobalTID(rtm, val) ()
#endif

#ifndef rtmGetMdlRefTriggerTID
# define rtmGetMdlRefTriggerTID(rtm)   ()
#endif

#ifndef rtmSetMdlRefTriggerTID
# define rtmSetMdlRefTriggerTID(rtm, val) ()
#endif

#ifndef rtmGetModelMappingInfo
# define rtmGetModelMappingInfo(rtm)   ((rtm)->SpecialInfo.mappingInfo)
#endif

#ifndef rtmSetModelMappingInfo
# define rtmSetModelMappingInfo(rtm, val) ((rtm)->SpecialInfo.mappingInfo = (val))
#endif

#ifndef rtmGetModelName
# define rtmGetModelName(rtm)          ((rtm)->modelName)
#endif

#ifndef rtmSetModelName
# define rtmSetModelName(rtm, val)     ((rtm)->modelName = (val))
#endif

#ifndef rtmGetNonInlinedSFcns
# define rtmGetNonInlinedSFcns(rtm)    ((rtm)->NonInlinedSFcns)
#endif

#ifndef rtmSetNonInlinedSFcns
# define rtmSetNonInlinedSFcns(rtm, val) ((rtm)->NonInlinedSFcns = (val))
#endif

#ifndef rtmGetNumBlockIO
# define rtmGetNumBlockIO(rtm)         ((rtm)->Sizes.numBlockIO)
#endif

#ifndef rtmSetNumBlockIO
# define rtmSetNumBlockIO(rtm, val)    ((rtm)->Sizes.numBlockIO = (val))
#endif

#ifndef rtmGetNumBlockParams
# define rtmGetNumBlockParams(rtm)     ((rtm)->Sizes.numBlockPrms)
#endif

#ifndef rtmSetNumBlockParams
# define rtmSetNumBlockParams(rtm, val) ((rtm)->Sizes.numBlockPrms = (val))
#endif

#ifndef rtmGetNumBlocks
# define rtmGetNumBlocks(rtm)          ((rtm)->Sizes.numBlocks)
#endif

#ifndef rtmSetNumBlocks
# define rtmSetNumBlocks(rtm, val)     ((rtm)->Sizes.numBlocks = (val))
#endif

#ifndef rtmGetNumContStates
# define rtmGetNumContStates(rtm)      ((rtm)->Sizes.numContStates)
#endif

#ifndef rtmSetNumContStates
# define rtmSetNumContStates(rtm, val) ((rtm)->Sizes.numContStates = (val))
#endif

#ifndef rtmGetNumDWork
# define rtmGetNumDWork(rtm)           ((rtm)->Sizes.numDwork)
#endif

#ifndef rtmSetNumDWork
# define rtmSetNumDWork(rtm, val)      ((rtm)->Sizes.numDwork = (val))
#endif

#ifndef rtmGetNumInputPorts
# define rtmGetNumInputPorts(rtm)      ((rtm)->Sizes.numIports)
#endif

#ifndef rtmSetNumInputPorts
# define rtmSetNumInputPorts(rtm, val) ((rtm)->Sizes.numIports = (val))
#endif

#ifndef rtmGetNumNonSampledZCs
# define rtmGetNumNonSampledZCs(rtm)   ((rtm)->Sizes.numNonSampZCs)
#endif

#ifndef rtmSetNumNonSampledZCs
# define rtmSetNumNonSampledZCs(rtm, val) ((rtm)->Sizes.numNonSampZCs = (val))
#endif

#ifndef rtmGetNumOutputPorts
# define rtmGetNumOutputPorts(rtm)     ((rtm)->Sizes.numOports)
#endif

#ifndef rtmSetNumOutputPorts
# define rtmSetNumOutputPorts(rtm, val) ((rtm)->Sizes.numOports = (val))
#endif

#ifndef rtmGetNumPeriodicContStates
# define rtmGetNumPeriodicContStates(rtm) ((rtm)->Sizes.numPeriodicContStates)
#endif

#ifndef rtmSetNumPeriodicContStates
# define rtmSetNumPeriodicContStates(rtm, val) ((rtm)->Sizes.numPeriodicContStates = (val))
#endif

#ifndef rtmGetNumSFcnParams
# define rtmGetNumSFcnParams(rtm)      ((rtm)->Sizes.numSFcnPrms)
#endif

#ifndef rtmSetNumSFcnParams
# define rtmSetNumSFcnParams(rtm, val) ((rtm)->Sizes.numSFcnPrms = (val))
#endif

#ifndef rtmGetNumSFunctions
# define rtmGetNumSFunctions(rtm)      ((rtm)->Sizes.numSFcns)
#endif

#ifndef rtmSetNumSFunctions
# define rtmSetNumSFunctions(rtm, val) ((rtm)->Sizes.numSFcns = (val))
#endif

#ifndef rtmGetNumSampleTimes
# define rtmGetNumSampleTimes(rtm)     ((rtm)->Sizes.numSampTimes)
#endif

#ifndef rtmSetNumSampleTimes
# define rtmSetNumSampleTimes(rtm, val) ((rtm)->Sizes.numSampTimes = (val))
#endif

#ifndef rtmGetNumU
# define rtmGetNumU(rtm)               ((rtm)->Sizes.numU)
#endif

#ifndef rtmSetNumU
# define rtmSetNumU(rtm, val)          ((rtm)->Sizes.numU = (val))
#endif

#ifndef rtmGetNumY
# define rtmGetNumY(rtm)               ((rtm)->Sizes.numY)
#endif

#ifndef rtmSetNumY
# define rtmSetNumY(rtm, val)          ((rtm)->Sizes.numY = (val))
#endif

#ifndef rtmGetOdeDELTA
# define rtmGetOdeDELTA(rtm)           ((rtm)->odeDELTA)
#endif

#ifndef rtmSetOdeDELTA
# define rtmSetOdeDELTA(rtm, val)      ((rtm)->odeDELTA = (val))
#endif

#ifndef rtmGetOdeDFDX
# define rtmGetOdeDFDX(rtm)            ((rtm)->odeDFDX)
#endif

#ifndef rtmSetOdeDFDX
# define rtmSetOdeDFDX(rtm, val)       ((rtm)->odeDFDX = (val))
#endif

#ifndef rtmGetOdeE
# define rtmGetOdeE(rtm)               ((rtm)->odeE)
#endif

#ifndef rtmSetOdeE
# define rtmSetOdeE(rtm, val)          ((rtm)->odeE = (val))
#endif

#ifndef rtmGetOdeF0
# define rtmGetOdeF0(rtm)              ((rtm)->odeF0)
#endif

#ifndef rtmSetOdeF0
# define rtmSetOdeF0(rtm, val)         ((rtm)->odeF0 = (val))
#endif

#ifndef rtmGetOdeF1
# define rtmGetOdeF1(rtm)              ((rtm)->odeF1)
#endif

#ifndef rtmSetOdeF1
# define rtmSetOdeF1(rtm, val)         ((rtm)->odeF1 = (val))
#endif

#ifndef rtmGetOdeFAC
# define rtmGetOdeFAC(rtm)             ((rtm)->odeFAC)
#endif

#ifndef rtmSetOdeFAC
# define rtmSetOdeFAC(rtm, val)        ((rtm)->odeFAC = (val))
#endif

#ifndef rtmGetOdePIVOTS
# define rtmGetOdePIVOTS(rtm)          ((rtm)->odePIVOTS)
#endif

#ifndef rtmSetOdePIVOTS
# define rtmSetOdePIVOTS(rtm, val)     ((rtm)->odePIVOTS = (val))
#endif

#ifndef rtmGetOdeW
# define rtmGetOdeW(rtm)               ((rtm)->odeW)
#endif

#ifndef rtmSetOdeW
# define rtmSetOdeW(rtm, val)          ((rtm)->odeW = (val))
#endif

#ifndef rtmGetOdeX0
# define rtmGetOdeX0(rtm)              ((rtm)->odeX0)
#endif

#ifndef rtmSetOdeX0
# define rtmSetOdeX0(rtm, val)         ((rtm)->odeX0 = (val))
#endif

#ifndef rtmGetOdeX1START
# define rtmGetOdeX1START(rtm)         ((rtm)->odeX1START)
#endif

#ifndef rtmSetOdeX1START
# define rtmSetOdeX1START(rtm, val)    ((rtm)->odeX1START = (val))
#endif

#ifndef rtmGetOdeXTMP
# define rtmGetOdeXTMP(rtm)            ((rtm)->odeXTMP)
#endif

#ifndef rtmSetOdeXTMP
# define rtmSetOdeXTMP(rtm, val)       ((rtm)->odeXTMP = (val))
#endif

#ifndef rtmGetOdeZTMP
# define rtmGetOdeZTMP(rtm)            ((rtm)->odeZTMP)
#endif

#ifndef rtmSetOdeZTMP
# define rtmSetOdeZTMP(rtm, val)       ((rtm)->odeZTMP = (val))
#endif

#ifndef rtmGetOffsetTimeArray
# define rtmGetOffsetTimeArray(rtm)    ((rtm)->Timing.offsetTimesArray)
#endif

#ifndef rtmSetOffsetTimeArray
# define rtmSetOffsetTimeArray(rtm, val) ((rtm)->Timing.offsetTimesArray = (val))
#endif

#ifndef rtmGetOffsetTimePtr
# define rtmGetOffsetTimePtr(rtm)      ((rtm)->Timing.offsetTimes)
#endif

#ifndef rtmSetOffsetTimePtr
# define rtmSetOffsetTimePtr(rtm, val) ((rtm)->Timing.offsetTimes = (val))
#endif

#ifndef rtmGetOptions
# define rtmGetOptions(rtm)            ((rtm)->Sizes.options)
#endif

#ifndef rtmSetOptions
# define rtmSetOptions(rtm, val)       ((rtm)->Sizes.options = (val))
#endif

#ifndef rtmGetParamIsMalloced
# define rtmGetParamIsMalloced(rtm)    ()
#endif

#ifndef rtmSetParamIsMalloced
# define rtmSetParamIsMalloced(rtm, val) ()
#endif

#ifndef rtmGetPath
# define rtmGetPath(rtm)               ((rtm)->path)
#endif

#ifndef rtmSetPath
# define rtmSetPath(rtm, val)          ((rtm)->path = (val))
#endif

#ifndef rtmGetPerTaskSampleHits
# define rtmGetPerTaskSampleHits(rtm)  ()
#endif

#ifndef rtmSetPerTaskSampleHits
# define rtmSetPerTaskSampleHits(rtm, val) ()
#endif

#ifndef rtmGetPerTaskSampleHitsArray
# define rtmGetPerTaskSampleHitsArray(rtm) ((rtm)->Timing.perTaskSampleHitsArray)
#endif

#ifndef rtmSetPerTaskSampleHitsArray
# define rtmSetPerTaskSampleHitsArray(rtm, val) ((rtm)->Timing.perTaskSampleHitsArray = (val))
#endif

#ifndef rtmGetPerTaskSampleHitsPtr
# define rtmGetPerTaskSampleHitsPtr(rtm) ((rtm)->Timing.perTaskSampleHits)
#endif

#ifndef rtmSetPerTaskSampleHitsPtr
# define rtmSetPerTaskSampleHitsPtr(rtm, val) ((rtm)->Timing.perTaskSampleHits = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
# define rtmGetPeriodicContStateIndices(rtm) ((rtm)->periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
# define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
# define rtmGetPeriodicContStateRanges(rtm) ((rtm)->periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
# define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->periodicContStateRanges = (val))
#endif

#ifndef rtmGetPrevZCSigState
# define rtmGetPrevZCSigState(rtm)     ((rtm)->prevZCSigState)
#endif

#ifndef rtmSetPrevZCSigState
# define rtmSetPrevZCSigState(rtm, val) ((rtm)->prevZCSigState = (val))
#endif

#ifndef rtmGetRTWExtModeInfo
# define rtmGetRTWExtModeInfo(rtm)     ((rtm)->extModeInfo)
#endif

#ifndef rtmSetRTWExtModeInfo
# define rtmSetRTWExtModeInfo(rtm, val) ((rtm)->extModeInfo = (val))
#endif

#ifndef rtmGetRTWGeneratedSFcn
# define rtmGetRTWGeneratedSFcn(rtm)   ((rtm)->Sizes.rtwGenSfcn)
#endif

#ifndef rtmSetRTWGeneratedSFcn
# define rtmSetRTWGeneratedSFcn(rtm, val) ((rtm)->Sizes.rtwGenSfcn = (val))
#endif

#ifndef rtmGetRTWLogInfo
# define rtmGetRTWLogInfo(rtm)         ((rtm)->rtwLogInfo)
#endif

#ifndef rtmSetRTWLogInfo
# define rtmSetRTWLogInfo(rtm, val)    ((rtm)->rtwLogInfo = (val))
#endif

#ifndef rtmGetRTWRTModelMethodsInfo
# define rtmGetRTWRTModelMethodsInfo(rtm) ()
#endif

#ifndef rtmSetRTWRTModelMethodsInfo
# define rtmSetRTWRTModelMethodsInfo(rtm, val) ()
#endif

#ifndef rtmGetRTWSfcnInfo
# define rtmGetRTWSfcnInfo(rtm)        ((rtm)->sfcnInfo)
#endif

#ifndef rtmSetRTWSfcnInfo
# define rtmSetRTWSfcnInfo(rtm, val)   ((rtm)->sfcnInfo = (val))
#endif

#ifndef rtmGetRTWSolverInfo
# define rtmGetRTWSolverInfo(rtm)      ((rtm)->solverInfo)
#endif

#ifndef rtmSetRTWSolverInfo
# define rtmSetRTWSolverInfo(rtm, val) ((rtm)->solverInfo = (val))
#endif

#ifndef rtmGetRTWSolverInfoPtr
# define rtmGetRTWSolverInfoPtr(rtm)   ((rtm)->solverInfoPtr)
#endif

#ifndef rtmSetRTWSolverInfoPtr
# define rtmSetRTWSolverInfoPtr(rtm, val) ((rtm)->solverInfoPtr = (val))
#endif

#ifndef rtmGetReservedForXPC
# define rtmGetReservedForXPC(rtm)     ((rtm)->SpecialInfo.xpcData)
#endif

#ifndef rtmSetReservedForXPC
# define rtmSetReservedForXPC(rtm, val) ((rtm)->SpecialInfo.xpcData = (val))
#endif

#ifndef rtmGetRootDWork
# define rtmGetRootDWork(rtm)          ((rtm)->dwork)
#endif

#ifndef rtmSetRootDWork
# define rtmSetRootDWork(rtm, val)     ((rtm)->dwork = (val))
#endif

#ifndef rtmGetSFunctions
# define rtmGetSFunctions(rtm)         ((rtm)->childSfunctions)
#endif

#ifndef rtmSetSFunctions
# define rtmSetSFunctions(rtm, val)    ((rtm)->childSfunctions = (val))
#endif

#ifndef rtmGetSampleHitArray
# define rtmGetSampleHitArray(rtm)     ((rtm)->Timing.sampleHitArray)
#endif

#ifndef rtmSetSampleHitArray
# define rtmSetSampleHitArray(rtm, val) ((rtm)->Timing.sampleHitArray = (val))
#endif

#ifndef rtmGetSampleHitPtr
# define rtmGetSampleHitPtr(rtm)       ((rtm)->Timing.sampleHits)
#endif

#ifndef rtmSetSampleHitPtr
# define rtmSetSampleHitPtr(rtm, val)  ((rtm)->Timing.sampleHits = (val))
#endif

#ifndef rtmGetSampleTimeArray
# define rtmGetSampleTimeArray(rtm)    ((rtm)->Timing.sampleTimesArray)
#endif

#ifndef rtmSetSampleTimeArray
# define rtmSetSampleTimeArray(rtm, val) ((rtm)->Timing.sampleTimesArray = (val))
#endif

#ifndef rtmGetSampleTimePtr
# define rtmGetSampleTimePtr(rtm)      ((rtm)->Timing.sampleTimes)
#endif

#ifndef rtmSetSampleTimePtr
# define rtmSetSampleTimePtr(rtm, val) ((rtm)->Timing.sampleTimes = (val))
#endif

#ifndef rtmGetSampleTimeTaskIDArray
# define rtmGetSampleTimeTaskIDArray(rtm) ((rtm)->Timing.sampleTimeTaskIDArray)
#endif

#ifndef rtmSetSampleTimeTaskIDArray
# define rtmSetSampleTimeTaskIDArray(rtm, val) ((rtm)->Timing.sampleTimeTaskIDArray = (val))
#endif

#ifndef rtmGetSampleTimeTaskIDPtr
# define rtmGetSampleTimeTaskIDPtr(rtm) ((rtm)->Timing.sampleTimeTaskIDPtr)
#endif

#ifndef rtmSetSampleTimeTaskIDPtr
# define rtmSetSampleTimeTaskIDPtr(rtm, val) ((rtm)->Timing.sampleTimeTaskIDPtr = (val))
#endif

#ifndef rtmGetSelf
# define rtmGetSelf(rtm)               ()
#endif

#ifndef rtmSetSelf
# define rtmSetSelf(rtm, val)          ()
#endif

#ifndef rtmGetSimMode
# define rtmGetSimMode(rtm)            ((rtm)->simMode)
#endif

#ifndef rtmSetSimMode
# define rtmSetSimMode(rtm, val)       ((rtm)->simMode = (val))
#endif

#ifndef rtmGetSimTimeStep
# define rtmGetSimTimeStep(rtm)        ((rtm)->Timing.simTimeStep)
#endif

#ifndef rtmSetSimTimeStep
# define rtmSetSimTimeStep(rtm, val)   ((rtm)->Timing.simTimeStep = (val))
#endif

#ifndef rtmGetStartTime
# define rtmGetStartTime(rtm)          ((rtm)->Timing.tStart)
#endif

#ifndef rtmSetStartTime
# define rtmSetStartTime(rtm, val)     ((rtm)->Timing.tStart = (val))
#endif

#ifndef rtmGetStepSize
# define rtmGetStepSize(rtm)           ((rtm)->Timing.stepSize)
#endif

#ifndef rtmSetStepSize
# define rtmSetStepSize(rtm, val)      ((rtm)->Timing.stepSize = (val))
#endif

#ifndef rtmGetStopRequestedFlag
# define rtmGetStopRequestedFlag(rtm)  ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequestedFlag
# define rtmSetStopRequestedFlag(rtm, val) ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetTaskCounters
# define rtmGetTaskCounters(rtm)       ((rtm)->Timing.TaskCounters)
#endif

#ifndef rtmSetTaskCounters
# define rtmSetTaskCounters(rtm, val)  ((rtm)->Timing.TaskCounters = (val))
#endif

#ifndef rtmGetTaskTimeArray
# define rtmGetTaskTimeArray(rtm)      ((rtm)->Timing.tArray)
#endif

#ifndef rtmSetTaskTimeArray
# define rtmSetTaskTimeArray(rtm, val) ((rtm)->Timing.tArray = (val))
#endif

#ifndef rtmGetTimePtr
# define rtmGetTimePtr(rtm)            ((rtm)->Timing.t)
#endif

#ifndef rtmSetTimePtr
# define rtmSetTimePtr(rtm, val)       ((rtm)->Timing.t = (val))
#endif

#ifndef rtmGetTimingData
# define rtmGetTimingData(rtm)         ((rtm)->Timing.timingData)
#endif

#ifndef rtmSetTimingData
# define rtmSetTimingData(rtm, val)    ((rtm)->Timing.timingData = (val))
#endif

#ifndef rtmGetU
# define rtmGetU(rtm)                  ((rtm)->inputs)
#endif

#ifndef rtmSetU
# define rtmSetU(rtm, val)             ((rtm)->inputs = (val))
#endif

#ifndef rtmGetVarNextHitTimesListPtr
# define rtmGetVarNextHitTimesListPtr(rtm) ((rtm)->Timing.varNextHitTimesList)
#endif

#ifndef rtmSetVarNextHitTimesListPtr
# define rtmSetVarNextHitTimesListPtr(rtm, val) ((rtm)->Timing.varNextHitTimesList = (val))
#endif

#ifndef rtmGetY
# define rtmGetY(rtm)                  ((rtm)->outputs)
#endif

#ifndef rtmSetY
# define rtmSetY(rtm, val)             ((rtm)->outputs = (val))
#endif

#ifndef rtmGetZCCacheNeedsReset
# define rtmGetZCCacheNeedsReset(rtm)  ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
# define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetZCSignalValues
# define rtmGetZCSignalValues(rtm)     ((rtm)->zcSignalValues)
#endif

#ifndef rtmSetZCSignalValues
# define rtmSetZCSignalValues(rtm, val) ((rtm)->zcSignalValues = (val))
#endif

#ifndef rtmGet_TimeOfLastOutput
# define rtmGet_TimeOfLastOutput(rtm)  ((rtm)->Timing.timeOfLastOutput)
#endif

#ifndef rtmSet_TimeOfLastOutput
# define rtmSet_TimeOfLastOutput(rtm, val) ((rtm)->Timing.timeOfLastOutput = (val))
#endif

#ifndef rtmGetdX
# define rtmGetdX(rtm)                 ((rtm)->derivs)
#endif

#ifndef rtmSetdX
# define rtmSetdX(rtm, val)            ((rtm)->derivs = (val))
#endif

#ifndef rtmGetChecksumVal
# define rtmGetChecksumVal(rtm, idx)   ((rtm)->Sizes.checksums[idx])
#endif

#ifndef rtmSetChecksumVal
# define rtmSetChecksumVal(rtm, idx, val) ((rtm)->Sizes.checksums[idx] = (val))
#endif

#ifndef rtmGetDWork
# define rtmGetDWork(rtm, idx)         ((rtm)->dwork[idx])
#endif

#ifndef rtmSetDWork
# define rtmSetDWork(rtm, idx, val)    ((rtm)->dwork[idx] = (val))
#endif

#ifndef rtmGetOffsetTime
# define rtmGetOffsetTime(rtm, idx)    ((rtm)->Timing.offsetTimes[idx])
#endif

#ifndef rtmSetOffsetTime
# define rtmSetOffsetTime(rtm, idx, val) ((rtm)->Timing.offsetTimes[idx] = (val))
#endif

#ifndef rtmGetSFunction
# define rtmGetSFunction(rtm, idx)     ((rtm)->childSfunctions[idx])
#endif

#ifndef rtmSetSFunction
# define rtmSetSFunction(rtm, idx, val) ((rtm)->childSfunctions[idx] = (val))
#endif

#ifndef rtmGetSampleTime
# define rtmGetSampleTime(rtm, idx)    ((rtm)->Timing.sampleTimes[idx])
#endif

#ifndef rtmSetSampleTime
# define rtmSetSampleTime(rtm, idx, val) ((rtm)->Timing.sampleTimes[idx] = (val))
#endif

#ifndef rtmGetSampleTimeTaskID
# define rtmGetSampleTimeTaskID(rtm, idx) ((rtm)->Timing.sampleTimeTaskIDPtr[idx])
#endif

#ifndef rtmSetSampleTimeTaskID
# define rtmSetSampleTimeTaskID(rtm, idx, val) ((rtm)->Timing.sampleTimeTaskIDPtr[idx] = (val))
#endif

#ifndef rtmGetVarNextHitTimeList
# define rtmGetVarNextHitTimeList(rtm, idx) ((rtm)->Timing.varNextHitTimesList[idx])
#endif

#ifndef rtmSetVarNextHitTimeList
# define rtmSetVarNextHitTimeList(rtm, idx, val) ((rtm)->Timing.varNextHitTimesList[idx] = (val))
#endif

#ifndef rtmIsContinuousTask
# define rtmIsContinuousTask(rtm, tid) ((tid) == 0)
#endif

#ifndef rtmGetErrorStatus
# define rtmGetErrorStatus(rtm)        ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
# define rtmSetErrorStatus(rtm, val)   ((rtm)->errorStatus = (val))
#endif

#ifndef rtmIsMajorTimeStep
# define rtmIsMajorTimeStep(rtm)       (((rtm)->Timing.simTimeStep) == MAJOR_TIME_STEP)
#endif

#ifndef rtmIsMinorTimeStep
# define rtmIsMinorTimeStep(rtm)       (((rtm)->Timing.simTimeStep) == MINOR_TIME_STEP)
#endif

#ifndef rtmIsSampleHit
# define rtmIsSampleHit(rtm, sti, tid) ((rtmIsMajorTimeStep((rtm)) && (rtm)->Timing.sampleHits[(rtm)->Timing.sampleTimeTaskIDPtr[sti]]))
#endif

#ifndef rtmGetStopRequested
# define rtmGetStopRequested(rtm)      ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
# define rtmSetStopRequested(rtm, val) ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
# define rtmGetStopRequestedPtr(rtm)   (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
# define rtmGetT(rtm)                  (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmSetT
# define rtmSetT(rtm, val)                                       /* Do Nothing */
#endif

#ifndef rtmGetTFinal
# define rtmGetTFinal(rtm)             ((rtm)->Timing.tFinal)
#endif

#ifndef rtmSetTFinal
# define rtmSetTFinal(rtm, val)        ((rtm)->Timing.tFinal = (val))
#endif

#ifndef rtmGetTPtr
# define rtmGetTPtr(rtm)               ((rtm)->Timing.t)
#endif

#ifndef rtmSetTPtr
# define rtmSetTPtr(rtm, val)          ((rtm)->Timing.t = (val))
#endif

#ifndef rtmGetTStart
# define rtmGetTStart(rtm)             ((rtm)->Timing.tStart)
#endif

#ifndef rtmSetTStart
# define rtmSetTStart(rtm, val)        ((rtm)->Timing.tStart = (val))
#endif

#ifndef rtmGetTaskTime
# define rtmGetTaskTime(rtm, sti)      (rtmGetTPtr((rtm))[(rtm)->Timing.sampleTimeTaskIDPtr[sti]])
#endif

#ifndef rtmSetTaskTime
# define rtmSetTaskTime(rtm, sti, val) (rtmGetTPtr((rtm))[sti] = (val))
#endif

#ifndef rtmGetTimeOfLastOutput
# define rtmGetTimeOfLastOutput(rtm)   ((rtm)->Timing.timeOfLastOutput)
#endif

#ifdef rtmGetRTWSolverInfo
#undef rtmGetRTWSolverInfo
#endif

#define rtmGetRTWSolverInfo(rtm)       &((rtm)->solverInfo)

/* Definition for use in the target main file */
#define mainV03_56_rtModel             RT_MODEL_mainV03_56_T

/* Block signals for system '<S47>/Distance into gust (y)' */
typedef struct {
  real_T DistanceintoGustxLimitedtogustlengthd;/* '<S51>/Distance into Gust (x) (Limited to gust length d) ' */
} B_Distanceintogusty_mainV03_56_T;

/* Block states (auto storage) for system '<S47>/Distance into gust (y)' */
typedef struct {
  boolean_T Distanceintogusty_MODE;    /* '<S47>/Distance into gust (y)' */
} DW_Distanceintogusty_mainV03_56_T;

/* Continuous states for system '<S47>/Distance into gust (y)' */
typedef struct {
  real_T DistanceintoGustxLimitedtogustlengthd_CSTATE;/* '<S51>/Distance into Gust (x) (Limited to gust length d) ' */
} X_Distanceintogusty_mainV03_56_T;

/* State derivatives for system '<S47>/Distance into gust (y)' */
typedef struct {
  real_T DistanceintoGustxLimitedtogustlengthd_CSTATE;/* '<S51>/Distance into Gust (x) (Limited to gust length d) ' */
} XDot_Distanceintogusty_mainV03_56_T;

/* State Disabled for system '<S47>/Distance into gust (y)' */
typedef struct {
  boolean_T DistanceintoGustxLimitedtogustlengthd_CSTATE;/* '<S51>/Distance into Gust (x) (Limited to gust length d) ' */
} XDis_Distanceintogusty_mainV03_56_T;

/* Zero-crossing (trigger) state for system '<S465>/J-K Flip-Flop' */
typedef struct {
  ZCSigState JKFlipFlop_Trig_ZCE;      /* '<S465>/J-K Flip-Flop' */
} ZCE_JKFlipFlop_mainV03_56_T;

/* Block signals (auto storage) */
typedef struct {
  plantDataBus plantData;              /* '<S4>/plantDataBus' */
  SensorsBus Sensors;                  /* '<S490>/Bus Creator2' */
  envDataBus envData;                  /* '<S2>/envDataBus' */
  analysisCasesBus analysisCasesBus_j; /* '<S242>/analysisCasesBus' */
  CommandBus BusCreator;               /* '<S5>/Bus Creator' */
  CommandBus BusCreator1;              /* '<S5>/Bus Creator1' */
  CommandBus ManualSwitch;             /* '<S5>/Manual Switch' */
  ThrottleBus Throttle;                /* '<S3>/Bus Creator1' */
  ThrottleBus ManualSwitch5;           /* '<S168>/Manual Switch5' */
  ThrottleBus BusCreator1_g;           /* '<S221>/Bus Creator1' */
  ThrottleBus ManualSwitch5_l;         /* '<S166>/Manual Switch5' */
  ThrottleBus ManualSwitch3;           /* '<S166>/Manual Switch3' */
  ThrottleBus BusCreator1_h;           /* '<S193>/Bus Creator1' */
  ThrottleBus ManualSwitch5_n;         /* '<S165>/Manual Switch5' */
  ThrottleBus ManualSwitch3_i;         /* '<S165>/Manual Switch3' */
  ThrottleBus ManualSwitch5_j;         /* '<S164>/Manual Switch5' */
  ThrottleBus Throttle_m;              /* '<S163>/Bus Creator1' */
  ActuatorsBus ActuatorsCmd;           /* '<S3>/Bus Creator' */
  ActuatorsBus actuators;              /* '<S168>/Bus Creator' */
  ActuatorsBus BusCreator1_l;          /* '<S168>/Bus Creator1' */
  ActuatorsBus ManualSwitch1;          /* '<S168>/Manual Switch1' */
  ActuatorsBus BusCreator2;            /* '<S166>/Bus Creator2' */
  ActuatorsBus actuators_h;            /* '<S165>/Bus Creator' */
  ActuatorsBus actuators_f;            /* '<S164>/Bus Creator' */
  ActuatorsBus BusCreator1_c;          /* '<S164>/Bus Creator1' */
  ActuatorsBus ManualSwitch1_n;        /* '<S164>/Manual Switch1' */
  ActuatorsBus Actuators;              /* '<S163>/Bus Creator' */
  LLABus LLA;                          /* '<S4>/LLABus' */
  EulerBus Euler;                      /* '<S4>/EulerBus' */
  OmegaBus Omega_body;                 /* '<S4>/OmegaBus' */
  betaBus betaBus_p;                   /* '<S242>/Bus Creator1' */
  betaBus betaBus_g;                   /* '<S242>/Bus Creator2' */
  alphaBus alphaBus_g;                 /* '<S242>/Bus Creator24' */
  alphaBus alphaBus_a;                 /* '<S242>/Bus Creator17' */
  altBus altBus_n;                     /* '<S242>/Bus Creator25' */
  altBus altBus_b;                     /* '<S242>/Bus Creator18' */
  deltaeBus deltaeBus_h;               /* '<S242>/Bus Creator27' */
  deltaeBus deltaeBus_f;               /* '<S242>/Bus Creator19' */
  deltarBus deltarBus_a;               /* '<S242>/Bus Creator28' */
  deltarBus deltarBus_b;               /* '<S242>/Bus Creator20' */
  deltafrBus deltafrBus_d;             /* '<S242>/Bus Creator29' */
  deltafrBus deltafrBus_l;             /* '<S242>/Bus Creator21' */
  deltaflBus deltaflBus_i;             /* '<S242>/Bus Creator30' */
  deltaflBus deltaflBus_k;             /* '<S242>/Bus Creator22' */
  xcgBus xcgBus_o;                     /* '<S242>/Bus Creator26' */
  xcgBus xcgBus_p;                     /* '<S242>/Bus Creator23' */
  real_T xeyeze[3];                    /* '<S243>/xe,ye,ze' */
  real_T SinCos_o1;                    /* '<S313>/SinCos' */
  real_T SinCos_o2;                    /* '<S313>/SinCos' */
  real_T xcos;                         /* '<S313>/x*cos' */
  real_T ysin;                         /* '<S313>/y*sin' */
  real_T Sum;                          /* '<S313>/Sum' */
  real_T Switch;                       /* '<S321>/Switch' */
  real_T TrigonometricFunction1;       /* '<S328>/Trigonometric Function1' */
  real_T radlat;                       /* '<S313>/rad lat' */
  real_T TrigonometricFunction2;       /* '<S328>/Trigonometric Function2' */
  real_T xsin;                         /* '<S313>/x*sin' */
  real_T ycos;                         /* '<S313>/y*cos' */
  real_T Sum1;                         /* '<S313>/Sum1' */
  real_T radlong;                      /* '<S313>/rad long ' */
  real_T Switch_a;                     /* '<S322>/Switch' */
  real_T Sum_g[2];                     /* '<S244>/Sum' */
  real_T Abs;                          /* '<S318>/Abs' */
  real_T Switch_f;                     /* '<S318>/Switch' */
  real_T Abs1;                         /* '<S315>/Abs1' */
  real_T Switch1;                      /* '<S311>/Switch1' */
  real_T Sum_c;                        /* '<S311>/Sum' */
  real_T Abs_h;                        /* '<S316>/Abs' */
  real_T phithetapsi[3];               /* '<S291>/phi theta psi' */
  real_T TmpSignalConversionAtsincosInport1[3];
  real_T VectorConcatenate[9];         /* '<S301>/Vector Concatenate' */
  real_T Transpose[9];                 /* '<S243>/Transpose' */
  real_T ubvbwb[3];                    /* '<S243>/ub,vb,wb' */
  real_T Product[3];                   /* '<S298>/Product' */
  real_T pqr[3];                       /* '<S243>/p,q,r ' */
  real_T TransferFcn;                  /* '<S250>/Transfer Fcn' */
  real_T TransferFcn1;                 /* '<S250>/Transfer Fcn1' */
  real_T TransferFcn2;                 /* '<S250>/Transfer Fcn2' */
  real_T dOmega_body[3];               /* '<S4>/TransferFunctions' */
  real_T TransferFcn3;                 /* '<S250>/Transfer Fcn3' */
  real_T TransferFcn4;                 /* '<S250>/Transfer Fcn4' */
  real_T TransferFcn5;                 /* '<S250>/Transfer Fcn5' */
  real_T Accel_body[3];                /* '<S4>/TransferFunctions' */
  real_T TransferFcn1_l;               /* '<S249>/Transfer Fcn1' */
  real_T TransferFcn4_i;               /* '<S249>/Transfer Fcn4' */
  real_T TransferFcn5_h;               /* '<S249>/Transfer Fcn5' */
  real_T Vaero[3];                     /* '<S4>/Sum2' */
  real_T alpha_j;                      /* '<S4>/IC7' */
  real_T Product_c;                    /* '<S335>/Product' */
  real_T Product1;                     /* '<S335>/Product1' */
  real_T Product2;                     /* '<S335>/Product2' */
  real_T Sum_o;                        /* '<S335>/Sum' */
  real_T Switch_h;                     /* '<S245>/Switch' */
  real_T beta_l;                       /* '<S4>/IC8' */
  real_T Derivative;                   /* '<S4>/Derivative' */
  real_T alpha_dot;                    /* '<S4>/Rate Limiter' */
  real_T Derivative1;                  /* '<S4>/Derivative1' */
  real_T beta_dot;                     /* '<S4>/Rate Limiter1' */
  real_T Sum3[3];                      /* '<S340>/Sum3' */
  real_T usedCapacityAs;               /* '<S247>/Integrator' */
  real_T usedCapacityAh;               /* '<S247>/Divide' */
  real_T nominalCapacityAh;            /* '<S247>/Product2' */
  real_T remainingCapacityAh;          /* '<S247>/Sum6' */
  real_T alpha_f;
  real_T beta_n;
  real_T Altitude_m;
  real_T X_ned[3];
  real_T phi;
  real_T theta;
  real_T psi;
  real_T Latitude_deg;
  real_T Longitude_deg;
  real_T Merge[4];                     /* '<S553>/Merge' */
  real_T V_body[3];
  real_T V_ned[3];
  real_T SFunction_o1;                 /* '<S12>/S-Function' */
  real_T SFunction_o2;                 /* '<S12>/S-Function' */
  real_T SFunction_o3;                 /* '<S12>/S-Function' */
  real_T SFunction_o4;                 /* '<S12>/S-Function' */
  real_T TmpSignalConversionAtSelectorInport1[3];
  real_T Selector2;                    /* '<S14>/Selector2' */
  real_T Selector1;                    /* '<S14>/Selector1' */
  real_T Selector;                     /* '<S14>/Selector' */
  real_T kclock;                       /* '<S2>/kclock' */
  real_T Sum1_j;                       /* '<S2>/Sum1' */
  real_T MatrixConcatenate[3];         /* '<S14>/Matrix Concatenate' */
  real_T GravityinEarthAxes[3];        /* '<S2>/Gravity in Earth Axes' */
  real_T MathFunction[3];              /* '<S2>/Math Function' */
  real_T Sum_p;                        /* '<S2>/Sum' */
  real_T r_A_body_1[3];                /* '<S13>/Sum2' */
  real_T ixj;                          /* '<S31>/i x j' */
  real_T jxk;                          /* '<S31>/j x k' */
  real_T kxi;                          /* '<S31>/k x i' */
  real_T ixk;                          /* '<S32>/i x k' */
  real_T jxi;                          /* '<S32>/j x i' */
  real_T kxj;                          /* '<S32>/k x j' */
  real_T Sum_e[3];                     /* '<S22>/Sum' */
  real_T Transpose_m[9];               /* '<S13>/Transpose' */
  real_T r_A_earth_1[3];               /* '<S13>/To Earth Axes2' */
  real_T r_A_earth_0[3];               /* '<S13>/Sum3' */
  real_T z;                            /* '<S28>/Sum5' */
  real_T Switch3;                      /* '<S28>/Switch3' */
  real_T K2;                           /* '<S28>/-K2' */
  real_T v_A_CG_0[3];                  /* '<S13>/To Earth Axes1' */
  real_T v_A_earth_0[3];               /* '<S13>/Sum4' */
  real_T D1;                           /* '<S28>/-D1' */
  real_T Sum4;                         /* '<S28>/Sum4' */
  real_T Product1_a;                   /* '<S28>/Product1' */
  real_T TmpSignalConversionAtToBodyAxes3Inport2[3];/* '<S13>/Reaction' */
  real_T F_A_body[3];                  /* '<S13>/To Body Axes3' */
  real_T ixj_j;                        /* '<S33>/i x j' */
  real_T jxk_d;                        /* '<S33>/j x k' */
  real_T kxi_c;                        /* '<S33>/k x i' */
  real_T ixk_n;                        /* '<S34>/i x k' */
  real_T jxi_g;                        /* '<S34>/j x i' */
  real_T kxj_m;                        /* '<S34>/k x j' */
  real_T Sum_d[3];                     /* '<S23>/Sum' */
  real_T r_A_body_1_g[3];              /* '<S13>/Sum6' */
  real_T ixj_i;                        /* '<S35>/i x j' */
  real_T jxk_p;                        /* '<S35>/j x k' */
  real_T kxi_i;                        /* '<S35>/k x i' */
  real_T ixk_o;                        /* '<S36>/i x k' */
  real_T jxi_f;                        /* '<S36>/j x i' */
  real_T kxj_i;                        /* '<S36>/k x j' */
  real_T Sum_k[3];                     /* '<S24>/Sum' */
  real_T r_A_earth_1_b[3];             /* '<S13>/To Earth Axes4' */
  real_T r_A_earth_0_n[3];             /* '<S13>/Sum7' */
  real_T z_l;                          /* '<S29>/Sum5' */
  real_T Switch3_k;                    /* '<S29>/Switch3' */
  real_T K2_n;                         /* '<S29>/-K2' */
  real_T v_A_CG_0_l[3];                /* '<S13>/To Earth Axes3' */
  real_T v_A_earth_0_b[3];             /* '<S13>/Sum8' */
  real_T D1_e;                         /* '<S29>/-D1' */
  real_T Sum4_j;                       /* '<S29>/Sum4' */
  real_T Product1_d;                   /* '<S29>/Product1' */
  real_T TmpSignalConversionAtToBodyAxes1Inport2[3];/* '<S13>/Reaction1' */
  real_T F_B_body[3];                  /* '<S13>/To Body Axes1' */
  real_T ixj_m;                        /* '<S37>/i x j' */
  real_T jxk_k;                        /* '<S37>/j x k' */
  real_T kxi_l;                        /* '<S37>/k x i' */
  real_T ixk_d;                        /* '<S38>/i x k' */
  real_T jxi_k;                        /* '<S38>/j x i' */
  real_T kxj_mr;                       /* '<S38>/k x j' */
  real_T Sum_a[3];                     /* '<S25>/Sum' */
  real_T r_A_body_1_e[3];              /* '<S13>/Sum9' */
  real_T ixj_b;                        /* '<S39>/i x j' */
  real_T jxk_c;                        /* '<S39>/j x k' */
  real_T kxi_g;                        /* '<S39>/k x i' */
  real_T ixk_c;                        /* '<S40>/i x k' */
  real_T jxi_j;                        /* '<S40>/j x i' */
  real_T kxj_a;                        /* '<S40>/k x j' */
  real_T Sum_cs[3];                    /* '<S26>/Sum' */
  real_T r_A_earth_1_i[3];             /* '<S13>/To Earth Axes6' */
  real_T r_A_earth_0_e[3];             /* '<S13>/Sum10' */
  real_T z_n;                          /* '<S30>/Sum5' */
  real_T Switch3_ks;                   /* '<S30>/Switch3' */
  real_T K2_l;                         /* '<S30>/-K2' */
  real_T v_A_CG_0_g[3];                /* '<S13>/To Earth Axes5' */
  real_T v_A_earth_0_h[3];             /* '<S13>/Sum11' */
  real_T D1_o;                         /* '<S30>/-D1' */
  real_T Sum4_f;                       /* '<S30>/Sum4' */
  real_T Product1_g;                   /* '<S30>/Product1' */
  real_T TmpSignalConversionAtToBodyAxes2Inport2[3];/* '<S13>/Reaction2' */
  real_T F_C_body[3];                  /* '<S13>/To Body Axes2' */
  real_T ixj_g;                        /* '<S41>/i x j' */
  real_T jxk_i;                        /* '<S41>/j x k' */
  real_T kxi_a;                        /* '<S41>/k x i' */
  real_T ixk_f;                        /* '<S42>/i x k' */
  real_T jxi_b;                        /* '<S42>/j x i' */
  real_T kxj_b;                        /* '<S42>/k x j' */
  real_T Sum_f[3];                     /* '<S27>/Sum' */
  real_T Sum_n;                        /* '<S104>/Sum' */
  real_T otime;                        /* '<S112>/otime' */
  real_T u80deg;                       /* '<S16>/+//- 180 deg' */
  real_T u0deg;                        /* '<S16>/+//- 90 deg' */
  real_T sincos_o1[2];                 /* '<S109>/sincos' */
  real_T sincos_o2[2];                 /* '<S109>/sincos' */
  real_T olon;                         /* '<S111>/olon' */
  real_T olat;                         /* '<S110>/olat' */
  real_T uto1000000m;                  /* '<S16>/0 to 1,000,000 m' */
  real_T Gain;                         /* '<S16>/Gain' */
  real_T oalt;                         /* '<S110>/oalt' */
  real_T Product_l;                    /* '<S109>/Product' */
  real_T Product1_m;                   /* '<S109>/Product1' */
  real_T aor;                          /* '<S104>/aor' */
  real_T by;                           /* '<S158>/Switch' */
  real_T Product1_aw;                  /* '<S157>/Product1' */
  real_T Product4;                     /* '<S157>/Product4' */
  real_T bx;                           /* '<S157>/Sum1' */
  real_T Product1_h;                   /* '<S159>/Product1' */
  real_T Product4_m;                   /* '<S159>/Product4' */
  real_T bz;                           /* '<S159>/Sum1' */
  real_T Product_k;                    /* '<S160>/Product' */
  real_T Product1_gj;                  /* '<S160>/Product1' */
  real_T Sum_pr;                       /* '<S160>/Sum' */
  real_T Product2_n;                   /* '<S160>/Product2' */
  real_T Sum1_p;                       /* '<S160>/Sum1' */
  real_T h1;                           /* '<S99>/h1' */
  real_T x1;                           /* '<S99>/x1' */
  real_T y1;                           /* '<S99>/y1' */
  real_T z1;                           /* '<S99>/z1' */
  real_T piGustlength[3];              /* '<S47>/pi//Gust length' */
  real_T Sum1_m[3];                    /* '<S47>/Sum1' */
  real_T Gustmagnitude20[3];           /* '<S47>/Gust magnitude//2.0' */
  real_T UnitConversion;               /* '<S54>/Unit Conversion' */
  real_T LimitFunction10ftto1000ft;    /* '<S91>/Limit Function 10ft to 1000ft' */
  real_T UnitConversion_j;             /* '<S93>/Unit Conversion' */
  real_T Lw[2];                        /* '<S62>/Lw' */
  real_T sigma_wg;                     /* '<S74>/sigma_wg ' */
  real_T PreLookUpIndexSearchaltitude_o2;/* '<S73>/PreLook-Up Index Search  (altitude)' */
  real_T PreLookUpIndexSearchprobofexceed_o2;/* '<S73>/PreLook-Up Index Search  (prob of exceed)' */
  real_T MediumHighAltitudeIntensity;  /* '<S73>/Medium//High Altitude Intensity' */
  real_T Output;                       /* '<S65>/Output' */
  real_T UnitConversion_l;             /* '<S58>/Unit Conversion' */
  real_T Airspeed;                     /* '<S15>/vt' */
  real_T Output_o[3];                  /* '<S66>/Output' */
  real_T LimitHeighth1000ft;           /* '<S74>/Limit Height h<1000ft' */
  real_T sigma_ugsigma_vg;             /* '<S74>/sigma_ug, sigma_vg' */
  real_T Lv[2];                        /* '<S62>/Lv' */
  real_T Merge_n[3];                   /* '<S78>/Merge' */
  real_T Merge_i[3];                   /* '<S86>/Merge' */
  real_T lnref_heightz0;               /* '<S49>/ln(ref_height//z0)' */
  real_T Windspeedatreferenceheight[3];/* '<S49>/Wind speed at reference height' */
  real_T axes[4];                      /* '<S5>/Joystick Input' */
  real_T POV;                          /* '<S5>/Joystick Input' */
  real_T ManualSwitch_m;               /* '<S3>/Manual Switch' */
  real_T Switch_e;                     /* '<S513>/Switch' */
  real_T Switch_n;                     /* '<S514>/Switch' */
  real_T Sum1_d[2];                    /* '<S502>/Sum1' */
  real_T Abs_k;                        /* '<S510>/Abs' */
  real_T Switch_ez;                    /* '<S510>/Switch' */
  real_T Abs1_j;                       /* '<S507>/Abs1' */
  real_T Switch_p;                     /* '<S507>/Switch' */
  real_T Switch1_g;                    /* '<S503>/Switch1' */
  real_T Sum_oa;                       /* '<S503>/Sum' */
  real_T Abs_d;                        /* '<S508>/Abs' */
  real_T Switch_aa;                    /* '<S508>/Switch' */
  real_T TrigonometricFunction1_k;     /* '<S520>/Trigonometric Function1' */
  real_T dNorth;                       /* '<S505>/dNorth' */
  real_T SinCos_o1_d;                  /* '<S505>/SinCos' */
  real_T SinCos_o2_p;                  /* '<S505>/SinCos' */
  real_T xcos_c;                       /* '<S505>/x*cos' */
  real_T TrigonometricFunction2_l;     /* '<S520>/Trigonometric Function2' */
  real_T dEast;                        /* '<S505>/dEast' */
  real_T ysin_l;                       /* '<S505>/y*sin' */
  real_T Px;                           /* '<S505>/Sum2' */
  real_T xsin_a;                       /* '<S505>/x*sin' */
  real_T ycos_b;                       /* '<S505>/y*cos' */
  real_T Py;                           /* '<S505>/Sum3' */
  real_T alt_p;                        /* '<S502>/Sum' */
  real_T X_nedMeas[3];                 /* '<S490>/GPS (Feedthough)' */
  real_T TransferFcnX;                 /* '<S538>/Transfer Fcn X' */
  real_T TransferFcnY;                 /* '<S538>/Transfer Fcn Y' */
  real_T TransferFcnZ;                 /* '<S538>/Transfer Fcn Z' */
  real_T ZeroOrderHold1[3];            /* '<S531>/Zero-Order Hold1' */
  real_T Product_b[3];                 /* '<S499>/Product' */
  real_T ZeroOrderHold2[3];            /* '<S531>/Zero-Order Hold2' */
  real_T ZeroOrderHold[3];             /* '<S531>/Zero-Order Hold' */
  real_T ZeroOrderHold4[3];            /* '<S531>/Zero-Order Hold4' */
  real_T Sum7[3];                      /* '<S531>/Sum7' */
  real_T Gain_d[3];                    /* '<S531>/Gain' */
  real_T jxk_b;                        /* '<S543>/j x k' */
  real_T kxi_f;                        /* '<S543>/k x i' */
  real_T ixj_e;                        /* '<S543>/i x j' */
  real_T kxj_l;                        /* '<S544>/k x j' */
  real_T ixk_o2;                       /* '<S544>/i x k' */
  real_T jxi_fh;                       /* '<S544>/j x i' */
  real_T Sum_b[3];                     /* '<S540>/Sum' */
  real_T jxk_cj;                       /* '<S541>/j x k' */
  real_T kxi_l1;                       /* '<S541>/k x i' */
  real_T ixj_k;                        /* '<S541>/i x j' */
  real_T kxj_me;                       /* '<S542>/k x j' */
  real_T ixk_fv;                       /* '<S542>/i x k' */
  real_T jxi_j3;                       /* '<S542>/j x i' */
  real_T Sum_h[3];                     /* '<S539>/Sum' */
  real_T ZeroOrderHold3[3];            /* '<S531>/Zero-Order Hold3' */
  real_T jxk_l;                        /* '<S545>/j x k' */
  real_T kxi_b;                        /* '<S545>/k x i' */
  real_T ixj_d;                        /* '<S545>/i x j' */
  real_T kxj_ly;                       /* '<S546>/k x j' */
  real_T ixk_h;                        /* '<S546>/i x k' */
  real_T jxi_ku;                       /* '<S546>/j x i' */
  real_T Sum_m[3];                     /* '<S536>/Sum' */
  real_T Sum_cg[3];                    /* '<S531>/Sum' */
  real_T Product_g[3];                 /* '<S531>/Product' */
  real_T Sum4_d[3];                    /* '<S531>/Sum4' */
  real_T Switch_o[3];                  /* '<S533>/Switch' */
  real_T Output_e[3];                  /* '<S534>/Output' */
  real_T Sum1_e[3];                    /* '<S531>/Sum1' */
  real_T TransferFcnX_k;               /* '<S550>/Transfer Fcn X' */
  real_T TransferFcnY_c;               /* '<S550>/Transfer Fcn Y' */
  real_T TransferFcnZ_b;               /* '<S550>/Transfer Fcn Z' */
  real_T ZeroOrderHold_m[3];           /* '<S532>/Zero-Order Hold' */
  real_T Product_a[3];                 /* '<S532>/Product' */
  real_T ZeroOrderHold1_f[3];          /* '<S532>/Zero-Order Hold1' */
  real_T Product1_aj[3];               /* '<S532>/Product1' */
  real_T Sum4_a[3];                    /* '<S532>/Sum4' */
  real_T Switch_l[3];                  /* '<S547>/Switch' */
  real_T Output_m[3];                  /* '<S548>/Output' */
  real_T Sum1_i[3];                    /* '<S532>/Sum1' */
  real_T Saturation_a[3];              /* '<S532>/Saturation' */
  real_T phithetapsi_m[3];             /* '<S525>/phi theta psi' */
  real_T TmpSignalConversionAtsincosInport1_p[3];
  real_T VectorConcatenate_h[9];       /* '<S529>/Vector Concatenate' */
  real_T Product7[4];                  /* '<S3>/Product7' */
  real_T Product5[4];                  /* '<S3>/Product5' */
  real_T Product3[4];                  /* '<S3>/Product3' */
  real_T Product1_m0[4];               /* '<S3>/Product1' */
  real_T Product9[4];                  /* '<S3>/Product9' */
  real_T ActuatorsCmd_b[4];            /* '<S3>/Sum1' */
  real_T Product6[5];                  /* '<S3>/Product6' */
  real_T Product4_o[5];                /* '<S3>/Product4' */
  real_T Product2_a[5];                /* '<S3>/Product2' */
  real_T Product_ki[5];                /* '<S3>/Product' */
  real_T Product8[5];                  /* '<S3>/Product8' */
  real_T Throttle_i[5];                /* '<S3>/Sum2' */
  real_T Sum_ba[3];                    /* '<S242>/Sum' */
  real_T VectorConcatenate_l[9];       /* '<S286>/Vector Concatenate' */
  real_T CD0;                          /* '<S253>/CD0' */
  real_T CDalpha;                      /* '<S253>/CDalpha' */
  real_T Product_d;                    /* '<S253>/Product' */
  real_T CDalpha_dot;                  /* '<S253>/CDalpha_dot' */
  real_T Product1_o;                   /* '<S253>/Product1' */
  real_T CDq;                          /* '<S253>/CDq' */
  real_T Product2_m;                   /* '<S253>/Product2' */
  real_T CDdeltae;                     /* '<S253>/CDdeltae' */
  real_T Product3_m;                   /* '<S253>/Product3' */
  real_T CDdeltafr;                    /* '<S253>/CDdeltafr' */
  real_T Product4_e;                   /* '<S253>/Product4' */
  real_T CDdeltafrl;                   /* '<S253>/CDdeltafrl' */
  real_T Product5_a;                   /* '<S253>/Product5' */
  real_T Add;                          /* '<S253>/Add' */
  real_T CYbeta;                       /* '<S255>/CYbeta' */
  real_T Product_gc;                   /* '<S255>/Product' */
  real_T CYp;                          /* '<S255>/CYp' */
  real_T Product1_p;                   /* '<S255>/Product1' */
  real_T CYr;                          /* '<S255>/CYr' */
  real_T Product2_j;                   /* '<S255>/Product2' */
  real_T CYdeltar;                     /* '<S255>/CYdeltar ' */
  real_T Product3_i;                   /* '<S255>/Product3' */
  real_T CYdeltafr;                    /* '<S255>/CYdeltafr' */
  real_T Product4_f;                   /* '<S255>/Product4' */
  real_T CYdeltafl;                    /* '<S255>/CYdeltafl' */
  real_T Product5_f;                   /* '<S255>/Product5' */
  real_T Add_g;                        /* '<S255>/Add' */
  real_T CL0;                          /* '<S254>/CL0' */
  real_T CLalpha;                      /* '<S254>/CLalpha' */
  real_T Product6_e;                   /* '<S254>/Product6' */
  real_T CLalpha_dot;                  /* '<S254>/CLalpha_dot' */
  real_T Product7_j;                   /* '<S254>/Product7' */
  real_T CLq;                          /* '<S254>/CLq' */
  real_T Product8_h;                   /* '<S254>/Product8' */
  real_T CLdeltae;                     /* '<S254>/CLdeltae' */
  real_T Product9_j;                   /* '<S254>/Product9' */
  real_T CLdeltafr;                    /* '<S254>/CLdeltafr' */
  real_T Product10;                    /* '<S254>/Product10' */
  real_T CLdeltafrl;                   /* '<S254>/CLdeltafrl' */
  real_T Product11;                    /* '<S254>/Product11' */
  real_T Add1;                         /* '<S254>/Add1' */
  real_T TmpSignalConversionAtProductInport2[3];
  real_T Product_i[3];                 /* '<S242>/Product' */
  real_T Product_f;                    /* '<S287>/Product' */
  real_T Product1_n;                   /* '<S287>/Product1' */
  real_T Product2_mn;                  /* '<S287>/Product2' */
  real_T Sum_c4;                       /* '<S287>/Sum' */
  real_T Product2_e;                   /* '<S260>/Product2' */
  real_T u2rhoV2;                      /* '<S260>/1//2rhoV^2' */
  real_T referencearea;                /* '<S252>/reference area' */
  real_T coefAdjust[3];                /* '<S252>/coefAdjust' */
  real_T Product_kd[3];                /* '<S252>/Product' */
  real_T VectorConcatenate_d[9];       /* '<S270>/Vector Concatenate' */
  real_T Sum_hf[3];                    /* '<S252>/Sum' */
  real_T Product_e[3];                 /* '<S263>/Product' */
  real_T ixj_p;                        /* '<S266>/i x j' */
  real_T jxk_do;                       /* '<S266>/j x k' */
  real_T kxi_gg;                       /* '<S266>/k x i' */
  real_T ixk_og;                       /* '<S267>/i x k' */
  real_T jxi_c;                        /* '<S267>/j x i' */
  real_T kxj_d;                        /* '<S267>/k x j' */
  real_T Sum_j[3];                     /* '<S262>/Sum' */
  real_T VectorConcatenate_d0[9];      /* '<S276>/Vector Concatenate' */
  real_T Transpose_h[9];               /* '<S264>/Transpose' */
  real_T Product_g0[3];                /* '<S264>/Product' */
  real_T VectorConcatenate_j[9];       /* '<S282>/Vector Concatenate' */
  real_T Transpose_e[9];               /* '<S265>/Transpose' */
  real_T Clbeta;                       /* '<S256>/Clbeta' */
  real_T Product6_h;                   /* '<S256>/Product6' */
  real_T Clp;                          /* '<S256>/Clp' */
  real_T Product7_e;                   /* '<S256>/Product7' */
  real_T Clr;                          /* '<S256>/Clr' */
  real_T Product8_j;                   /* '<S256>/Product8' */
  real_T Cldeltar;                     /* '<S256>/Cldeltar ' */
  real_T Product9_d;                   /* '<S256>/Product9' */
  real_T Cldeltafr;                    /* '<S256>/Cldeltafr' */
  real_T Product10_d;                  /* '<S256>/Product10' */
  real_T Cldeltafl;                    /* '<S256>/Cldeltafl' */
  real_T Product11_d;                  /* '<S256>/Product11' */
  real_T Add1_f;                       /* '<S256>/Add1' */
  real_T Cm0;                          /* '<S257>/Cm0' */
  real_T Cmalpha;                      /* '<S257>/Cmalpha' */
  real_T Product6_ep;                  /* '<S257>/Product6' */
  real_T Cmalpha_dot;                  /* '<S257>/Cmalpha_dot' */
  real_T Product7_k;                   /* '<S257>/Product7' */
  real_T Cmq;                          /* '<S257>/Cmq' */
  real_T Product8_p;                   /* '<S257>/Product8' */
  real_T Cmdeltae;                     /* '<S257>/Cmdeltae' */
  real_T Product9_p;                   /* '<S257>/Product9' */
  real_T Cmdeltafr;                    /* '<S257>/Cmdeltafr' */
  real_T Product10_j;                  /* '<S257>/Product10' */
  real_T Cmdeltafl;                    /* '<S257>/Cmdeltafl' */
  real_T Product11_a;                  /* '<S257>/Product11' */
  real_T Add2;                         /* '<S257>/Add2' */
  real_T Cnbeta;                       /* '<S258>/Cnbeta' */
  real_T Product6_m;                   /* '<S258>/Product6' */
  real_T Cnp;                          /* '<S258>/Cnp' */
  real_T Product7_f;                   /* '<S258>/Product7' */
  real_T Cnr;                          /* '<S258>/Cnr' */
  real_T Product8_f;                   /* '<S258>/Product8' */
  real_T Cndeltar;                     /* '<S258>/Cndeltar ' */
  real_T Product9_e;                   /* '<S258>/Product9' */
  real_T Cndeltafr;                    /* '<S258>/Cndeltafr' */
  real_T Product10_k;                  /* '<S258>/Product10' */
  real_T Cndeltafl;                    /* '<S258>/Cndeltafl' */
  real_T Product11_g;                  /* '<S258>/Product11' */
  real_T Add1_e;                       /* '<S258>/Add1' */
  real_T Product1_k[3];                /* '<S252>/Product1' */
  real_T Product3_ix[3];               /* '<S252>/Product3' */
  real_T Sum1_g[3];                    /* '<S252>/Sum1' */
  real_T Product_fr[3];                /* '<S265>/Product' */
  real_T phidot;                       /* '<S300>/phidot' */
  real_T psidot;                       /* '<S300>/psidot' */
  real_T thetadot;                     /* '<S300>/thetadot' */
  real_T Selector_j[9];                /* '<S292>/Selector' */
  real_T Product_fi[3];                /* '<S302>/Product' */
  real_T Selector1_d[9];               /* '<S292>/Selector1' */
  real_T Product_ca[3];                /* '<S303>/Product' */
  real_T arm[3];                       /* '<S349>/Sum1' */
  real_T jxk_g;                        /* '<S355>/j x k' */
  real_T Product_p;                    /* '<S352>/Product' */
  real_T Product1_kx;                  /* '<S352>/Product1' */
  real_T VoltagePrelookup_o2;          /* '<S353>/Voltage Prelookup' */
  real_T HeightPrelookup_o2;           /* '<S353>/Height Prelookup' */
  real_T Sum_a5;                       /* '<S349>/Sum' */
  real_T VflightPrelookup_o2;          /* '<S353>/Vflight Prelookup' */
  real_T Thrust;                       /* '<S353>/Thrust Interpolation Using Prelookup' */
  real_T kxi_m;                        /* '<S355>/k x i' */
  real_T ixj_n;                        /* '<S355>/i x j' */
  real_T kxj_g;                        /* '<S356>/k x j' */
  real_T ixk_a;                        /* '<S356>/i x k' */
  real_T jxi_bv;                       /* '<S356>/j x i' */
  real_T Sum_ac[3];                    /* '<S351>/Sum' */
  real_T RPM;                          /* '<S353>/RPM Interpolation Using Prelookup' */
  real_T Sum_kp;                       /* '<S354>/Sum' */
  real_T Torque;                       /* '<S354>/Switch' */
  real_T Sum2[3];                      /* '<S349>/Sum2' */
  real_T arm_h[3];                     /* '<S363>/Sum1' */
  real_T Product_dq;                   /* '<S366>/Product' */
  real_T Product1_c;                   /* '<S366>/Product1' */
  real_T VoltagePrelookup_o2_j;        /* '<S367>/Voltage Prelookup' */
  real_T HeightPrelookup_o2_k;         /* '<S367>/Height Prelookup' */
  real_T Sum_eo;                       /* '<S363>/Sum' */
  real_T VflightPrelookup_o2_d;        /* '<S367>/Vflight Prelookup' */
  real_T Thrust_d;                     /* '<S367>/Thrust Interpolation Using Prelookup' */
  real_T jxk_o;                        /* '<S369>/j x k' */
  real_T kxi_gm;                       /* '<S369>/k x i' */
  real_T ixj_c;                        /* '<S369>/i x j' */
  real_T kxj_gt;                       /* '<S370>/k x j' */
  real_T ixk_b;                        /* '<S370>/i x k' */
  real_T jxi_gj;                       /* '<S370>/j x i' */
  real_T Sum_h5[3];                    /* '<S365>/Sum' */
  real_T RPM_g;                        /* '<S367>/RPM Interpolation Using Prelookup' */
  real_T Sum_m4;                       /* '<S368>/Sum' */
  real_T Torque_h;                     /* '<S368>/Switch' */
  real_T Sum2_a[3];                    /* '<S363>/Sum2' */
  real_T arm_e[3];                     /* '<S377>/Sum1' */
  real_T Product_j;                    /* '<S380>/Product' */
  real_T Product1_e;                   /* '<S380>/Product1' */
  real_T VoltagePrelookup_o2_a;        /* '<S381>/Voltage Prelookup' */
  real_T HeightPrelookup_o2_j;         /* '<S381>/Height Prelookup' */
  real_T Sum_nq;                       /* '<S377>/Sum' */
  real_T VflightPrelookup_o2_f;        /* '<S381>/Vflight Prelookup' */
  real_T Thrust_l;                     /* '<S381>/Thrust Interpolation Using Prelookup' */
  real_T jxk_b2;                       /* '<S383>/j x k' */
  real_T kxi_go;                       /* '<S383>/k x i' */
  real_T ixj_dn;                       /* '<S383>/i x j' */
  real_T kxj_e;                        /* '<S384>/k x j' */
  real_T ixk_op;                       /* '<S384>/i x k' */
  real_T jxi_d;                        /* '<S384>/j x i' */
  real_T Sum_j1[3];                    /* '<S379>/Sum' */
  real_T RPM_f;                        /* '<S381>/RPM Interpolation Using Prelookup' */
  real_T Sum_cl;                       /* '<S382>/Sum' */
  real_T Torque_n;                     /* '<S382>/Switch' */
  real_T Sum2_h[3];                    /* '<S377>/Sum2' */
  real_T arm_p[3];                     /* '<S391>/Sum1' */
  real_T Product_ch;                   /* '<S394>/Product' */
  real_T Product1_kd;                  /* '<S394>/Product1' */
  real_T VoltagePrelookup_o2_k;        /* '<S395>/Voltage Prelookup' */
  real_T HeightPrelookup_o2_e;         /* '<S395>/Height Prelookup' */
  real_T Sum_fv;                       /* '<S391>/Sum' */
  real_T VflightPrelookup_o2_fy;       /* '<S395>/Vflight Prelookup' */
  real_T Thrust_e;                     /* '<S395>/Thrust Interpolation Using Prelookup' */
  real_T jxk_a;                        /* '<S397>/j x k' */
  real_T kxi_m5;                       /* '<S397>/k x i' */
  real_T ixj_n0;                       /* '<S397>/i x j' */
  real_T kxj_p;                        /* '<S398>/k x j' */
  real_T ixk_h1;                       /* '<S398>/i x k' */
  real_T jxi_m;                        /* '<S398>/j x i' */
  real_T Sum_l[3];                     /* '<S393>/Sum' */
  real_T RPM_i;                        /* '<S395>/RPM Interpolation Using Prelookup' */
  real_T Sum_pt;                       /* '<S396>/Sum' */
  real_T Torque_c;                     /* '<S396>/Switch' */
  real_T Sum2_m[3];                    /* '<S391>/Sum2' */
  real_T arm_a[3];                     /* '<S405>/Sum1' */
  real_T Product_h;                    /* '<S408>/Product' */
  real_T Product1_ef;                  /* '<S408>/Product1' */
  real_T VoltagePrelookup_o2_f;        /* '<S409>/Voltage Prelookup' */
  real_T HeightPrelookup_o2_o;         /* '<S409>/Height Prelookup' */
  real_T Sum_ar;                       /* '<S405>/Sum' */
  real_T VflightPrelookup_o2_m;        /* '<S409>/Vflight Prelookup' */
  real_T Thrust_m;                     /* '<S409>/Thrust Interpolation Using Prelookup' */
  real_T jxk_e;                        /* '<S411>/j x k' */
  real_T kxi_fl;                       /* '<S411>/k x i' */
  real_T ixj_eu;                       /* '<S411>/i x j' */
  real_T kxj_m5;                       /* '<S412>/k x j' */
  real_T ixk_nc;                       /* '<S412>/i x k' */
  real_T jxi_b4;                       /* '<S412>/j x i' */
  real_T Sum_ab[3];                    /* '<S407>/Sum' */
  real_T RPM_k;                        /* '<S409>/RPM Interpolation Using Prelookup' */
  real_T Sum_b0;                       /* '<S410>/Sum' */
  real_T Torque_p;                     /* '<S410>/Switch' */
  real_T Sum2_c[3];                    /* '<S405>/Sum2' */
  real_T total_moment[3];              /* '<S247>/Sum3' */
  real_T Sum1_p3[3];                   /* '<S4>/Sum1' */
  real_T jxk_cb;                       /* '<S305>/j x k' */
  real_T kxi_j;                        /* '<S305>/k x i' */
  real_T ixj_ey;                       /* '<S305>/i x j' */
  real_T kxj_d5;                       /* '<S306>/k x j' */
  real_T ixk_fm;                       /* '<S306>/i x k' */
  real_T jxi_o;                        /* '<S306>/j x i' */
  real_T Sum_bh[3];                    /* '<S304>/Sum' */
  real_T Sum2_mg[3];                   /* '<S292>/Sum2' */
  real_T Selector2_n[9];               /* '<S292>/Selector2' */
  real_T Product2_ex[3];               /* '<S292>/Product2' */
  real_T SumofElements[3];             /* '<S307>/Sum of Elements' */
  real_T total_thrust[3];              /* '<S247>/Sum2' */
  real_T Sum_ma[3];                    /* '<S4>/Sum' */
  real_T Sum_ny[3];                    /* '<S293>/Sum' */
  real_T Product_km[3];                /* '<S243>/Product' */
  real_T jxk_j;                        /* '<S309>/j x k' */
  real_T kxi_k;                        /* '<S309>/k x i' */
  real_T ixj_f;                        /* '<S309>/i x j' */
  real_T kxj_n;                        /* '<S310>/k x j' */
  real_T ixk_bi;                       /* '<S310>/i x k' */
  real_T jxi_dx;                       /* '<S310>/j x i' */
  real_T Sum_as[3];                    /* '<S294>/Sum' */
  real_T Sum2_af[3];                   /* '<S243>/Sum2' */
  real_T dOmega_body_k[3];             /* '<S4>/IC' */
  real_T Accel_body_a[3];              /* '<S4>/IC1' */
  real_T Current;                      /* '<S353>/Current Interpolation Using Prelookup' */
  real_T Current_d;                    /* '<S367>/Current Interpolation Using Prelookup' */
  real_T Current_b;                    /* '<S381>/Current Interpolation Using Prelookup' */
  real_T Current_m;                    /* '<S395>/Current Interpolation Using Prelookup' */
  real_T Current_h;                    /* '<S409>/Current Interpolation Using Prelookup' */
  real_T total_Im;                     /* '<S247>/Sum5' */
  real_T Sum_db;                       /* '<S464>/Sum' */
  real_T phidot_p;                     /* '<S528>/phidot' */
  real_T psidot_k;                     /* '<S528>/psidot' */
  real_T thetadot_i;                   /* '<S528>/thetadot' */
  real_T Subtract;                     /* '<S556>/Subtract' */
  real_T Gain1;                        /* '<S556>/Gain1' */
  real_T Sum_eh[3];                    /* '<S556>/Sum' */
  real_T Gain_a;                       /* '<S556>/Gain' */
  real_T Sum_bt[3];                    /* '<S558>/Sum' */
  real_T Gain_c[3];                    /* '<S558>/Gain' */
  real_T DeadZone[3];                  /* '<S559>/Dead Zone' */
  real_T Sign1[3];                     /* '<S559>/Sign1' */
  real_T SumofElements3;               /* '<S559>/Sum of Elements3' */
  real_T Switch_dl[3];                 /* '<S559>/Switch' */
  real_T MinMax[3];                    /* '<S558>/MinMax' */
  real_T Switch_dj[3];                 /* '<S558>/Switch' */
  real_T Product_o[3];                 /* '<S558>/Product' */
  real_T ProductofElements;            /* '<S559>/Product of Elements' */
  real_T Switch1_a[3];                 /* '<S559>/Switch1' */
  real_T Gain1_m[3];                   /* '<S559>/Gain1' */
  real_T Bias;                         /* '<S508>/Bias' */
  real_T Bias1;                        /* '<S508>/Bias1' */
  real_T Bias_a;                       /* '<S510>/Bias' */
  real_T Bias1_h;                      /* '<S510>/Bias1' */
  real_T Bias_o;                       /* '<S507>/Bias' */
  real_T Gain_j;                       /* '<S507>/Gain' */
  real_T Bias1_hf;                     /* '<S507>/Bias1' */
  real_T Sign1_h;                      /* '<S507>/Sign1' */
  real_T Divide1;                      /* '<S507>/Divide1' */
  real_T VflightPrelookup_o2_n;        /* '<S410>/Vflight Prelookup' */
  real_T RPMPrelookup_o2;              /* '<S410>/RPM Prelookup' */
  real_T Torque_b;                     /* '<S410>/Torque Interpolation Using Prelookup' */
  real_T VflightPrelookup_o2_k;        /* '<S396>/Vflight Prelookup' */
  real_T RPMPrelookup_o2_e;            /* '<S396>/RPM Prelookup' */
  real_T Torque_pl;                    /* '<S396>/Torque Interpolation Using Prelookup' */
  real_T VflightPrelookup_o2_a;        /* '<S382>/Vflight Prelookup' */
  real_T RPMPrelookup_o2_a;            /* '<S382>/RPM Prelookup' */
  real_T Torque_hw;                    /* '<S382>/Torque Interpolation Using Prelookup' */
  real_T VflightPrelookup_o2_j;        /* '<S368>/Vflight Prelookup' */
  real_T RPMPrelookup_o2_ao;           /* '<S368>/RPM Prelookup' */
  real_T Torque_a;                     /* '<S368>/Torque Interpolation Using Prelookup' */
  real_T VflightPrelookup_o2_dh;       /* '<S354>/Vflight Prelookup' */
  real_T RPMPrelookup_o2_n;            /* '<S354>/RPM Prelookup' */
  real_T Torque_az;                    /* '<S354>/Torque Interpolation Using Prelookup' */
  real_T Product_m;                    /* '<S245>/Product' */
  real_T Bias_e;                       /* '<S316>/Bias' */
  real_T Bias1_e;                      /* '<S316>/Bias1' */
  real_T Bias_b;                       /* '<S318>/Bias' */
  real_T Bias1_c;                      /* '<S318>/Bias1' */
  real_T Bias_l;                       /* '<S315>/Bias' */
  real_T Gain_o;                       /* '<S315>/Gain' */
  real_T Bias1_g;                      /* '<S315>/Bias1' */
  real_T Sign1_f;                      /* '<S315>/Sign1' */
  real_T Divide1_o;                    /* '<S315>/Divide1' */
  real_T Sum1_b;                       /* '<S168>/Sum1' */
  real_T Saturation_c;                 /* '<S168>/Saturation' */
  real_T ProportionalGain;             /* '<S233>/Proportional Gain' */
  real_T Integrator;                   /* '<S233>/Integrator' */
  real_T DerivativeGain;               /* '<S233>/Derivative Gain' */
  real_T Filter;                       /* '<S233>/Filter' */
  real_T SumD;                         /* '<S233>/SumD' */
  real_T FilterCoefficient;            /* '<S233>/Filter Coefficient' */
  real_T Sum_jd;                       /* '<S233>/Sum' */
  real_T Saturate;                     /* '<S233>/Saturate' */
  real_T Sum12;                        /* '<S168>/Sum12' */
  real_T ProportionalGain_f;           /* '<S235>/Proportional Gain' */
  real_T Integrator_l;                 /* '<S235>/Integrator' */
  real_T DerivativeGain_d;             /* '<S235>/Derivative Gain' */
  real_T Filter_h;                     /* '<S235>/Filter' */
  real_T SumD_f;                       /* '<S235>/SumD' */
  real_T FilterCoefficient_l;          /* '<S235>/Filter Coefficient' */
  real_T Sum_pn;                       /* '<S235>/Sum' */
  real_T PitchRateIn;                  /* '<S168>/Sum13' */
  real_T ProportionalGain_n;           /* '<S234>/Proportional Gain' */
  real_T Integrator_e;                 /* '<S234>/Integrator' */
  real_T DerivativeGain_f;             /* '<S234>/Derivative Gain' */
  real_T Filter_f;                     /* '<S234>/Filter' */
  real_T SumD_i;                       /* '<S234>/SumD' */
  real_T FilterCoefficient_c;          /* '<S234>/Filter Coefficient' */
  real_T Sum_g1;                       /* '<S234>/Sum' */
  real_T Saturate_i;                   /* '<S234>/Saturate' */
  real_T deltae_k;                     /* '<S168>/Gain4' */
  real_T Sum3_j;                       /* '<S168>/Sum3' */
  real_T ProportionalGain_l;           /* '<S237>/Proportional Gain' */
  real_T Integrator_f;                 /* '<S237>/Integrator' */
  real_T DerivativeGain_a;             /* '<S237>/Derivative Gain' */
  real_T Filter_d;                     /* '<S237>/Filter' */
  real_T SumD_ig;                      /* '<S237>/SumD' */
  real_T FilterCoefficient_k;          /* '<S237>/Filter Coefficient' */
  real_T Sum_jk;                       /* '<S237>/Sum' */
  real_T Sum2_i;                       /* '<S168>/Sum2' */
  real_T ProportionalGain_g;           /* '<S236>/Proportional Gain' */
  real_T Integrator_g;                 /* '<S236>/Integrator' */
  real_T DerivativeGain_h;             /* '<S236>/Derivative Gain' */
  real_T Filter_dv;                    /* '<S236>/Filter' */
  real_T SumD_c;                       /* '<S236>/SumD' */
  real_T FilterCoefficient_o;          /* '<S236>/Filter Coefficient' */
  real_T deltar_m;                     /* '<S168>/Gain5' */
  real_T Sum5_h;                       /* '<S168>/Sum5' */
  real_T ProportionalGain_o;           /* '<S239>/Proportional Gain' */
  real_T Integrator_lp;                /* '<S239>/Integrator' */
  real_T DerivativeGain_e;             /* '<S239>/Derivative Gain' */
  real_T Filter_fj;                    /* '<S239>/Filter' */
  real_T SumD_l;                       /* '<S239>/SumD' */
  real_T FilterCoefficient_km;         /* '<S239>/Filter Coefficient' */
  real_T Sum_ck;                       /* '<S239>/Sum' */
  real_T Sum4_h;                       /* '<S168>/Sum4' */
  real_T ProportionalGain_a;           /* '<S238>/Proportional Gain' */
  real_T Integrator_lx;                /* '<S238>/Integrator' */
  real_T DerivativeGain_j;             /* '<S238>/Derivative Gain' */
  real_T Filter_dd;                    /* '<S238>/Filter' */
  real_T SumD_h;                       /* '<S238>/SumD' */
  real_T FilterCoefficient_m;          /* '<S238>/Filter Coefficient' */
  real_T Sum_f0;                       /* '<S238>/Sum' */
  real_T deltafr_j;                    /* '<S168>/Gain6' */
  real_T deltafl_o;                    /* '<S168>/Gain7' */
  real_T IntegralGain;                 /* '<S233>/Integral Gain' */
  real_T IntegralGain_d;               /* '<S234>/Integral Gain' */
  real_T IntegralGain_i;               /* '<S235>/Integral Gain' */
  real_T IntegralGain_f;               /* '<S236>/Integral Gain' */
  real_T IntegralGain_a;               /* '<S237>/Integral Gain' */
  real_T IntegralGain_j;               /* '<S238>/Integral Gain' */
  real_T IntegralGain_b;               /* '<S239>/Integral Gain' */
  real_T Sum_jn;                       /* '<S236>/Sum' */
  real_T Sum6;                         /* '<S168>/Sum6' */
  real_T Gain8;                        /* '<S168>/Gain8' */
  real_T Sum_hw;                       /* '<S168>/Sum' */
  real_T TmpSignalConversionAtsincosInport1_f[3];
  real_T VectorConcatenate_e[9];       /* '<S227>/Vector Concatenate' */
  real_T ToBodyAxes3[3];               /* '<S166>/To Body Axes3' */
  real_T ManualSwitch10;               /* '<S166>/Manual Switch10' */
  real_T Sum18;                        /* '<S166>/Sum18' */
  real_T Saturation11;                 /* '<S166>/Saturation11' */
  real_T ProportionalGain_k;           /* '<S201>/Proportional Gain' */
  real_T Integrator_gn;                /* '<S201>/Integrator' */
  real_T DerivativeGain_da;            /* '<S201>/Derivative Gain' */
  real_T Filter_fi;                    /* '<S201>/Filter' */
  real_T SumD_a;                       /* '<S201>/SumD' */
  real_T FilterCoefficient_j;          /* '<S201>/Filter Coefficient' */
  real_T Sum_pru;                      /* '<S201>/Sum' */
  real_T Saturate_m;                   /* '<S201>/Saturate' */
  real_T Sum19;                        /* '<S166>/Sum19' */
  real_T ProportionalGain_kl;          /* '<S207>/Proportional Gain' */
  real_T Integrator_i;                 /* '<S207>/Integrator' */
  real_T DerivativeGain_g;             /* '<S207>/Derivative Gain' */
  real_T Filter_o;                     /* '<S207>/Filter' */
  real_T SumD_ie;                      /* '<S207>/SumD' */
  real_T FilterCoefficient_f;          /* '<S207>/Filter Coefficient' */
  real_T Sum_a2;                       /* '<S207>/Sum' */
  real_T PitchRateIn_j;                /* '<S166>/Sum20' */
  real_T ProportionalGain_nu;          /* '<S206>/Proportional Gain' */
  real_T Integrator_n;                 /* '<S206>/Integrator' */
  real_T DerivativeGain_i;             /* '<S206>/Derivative Gain' */
  real_T Filter_e;                     /* '<S206>/Filter' */
  real_T SumD_g;                       /* '<S206>/SumD' */
  real_T FilterCoefficient_ln;         /* '<S206>/Filter Coefficient' */
  real_T Sum_ph;                       /* '<S206>/Sum' */
  real_T Saturate_k;                   /* '<S206>/Saturate' */
  real_T deltae_p;                     /* '<S166>/Gain4' */
  real_T deltae_l;                     /* '<S166>/Manual Switch2' */
  real_T Sum22;                        /* '<S166>/Sum22' */
  real_T ProportionalGain_gi;          /* '<S209>/Proportional Gain' */
  real_T Integrator_en;                /* '<S209>/Integrator' */
  real_T DerivativeGain_eh;            /* '<S209>/Derivative Gain' */
  real_T Filter_l;                     /* '<S209>/Filter' */
  real_T SumD_ch;                      /* '<S209>/SumD' */
  real_T FilterCoefficient_mb;         /* '<S209>/Filter Coefficient' */
  real_T Sum_mw;                       /* '<S209>/Sum' */
  real_T Sum21;                        /* '<S166>/Sum21' */
  real_T ProportionalGain_c;           /* '<S208>/Proportional Gain' */
  real_T Integrator_d;                 /* '<S208>/Integrator' */
  real_T DerivativeGain_n;             /* '<S208>/Derivative Gain' */
  real_T Filter_d2;                    /* '<S208>/Filter' */
  real_T SumD_p;                       /* '<S208>/SumD' */
  real_T FilterCoefficient_n;          /* '<S208>/Filter Coefficient' */
  real_T deltar_p;                     /* '<S166>/Gain5' */
  real_T deltar_j;                     /* '<S166>/Manual Switch4' */
  real_T Sum24;                        /* '<S166>/Sum24' */
  real_T ProportionalGain_b;           /* '<S211>/Proportional Gain' */
  real_T Integrator_b;                 /* '<S211>/Integrator' */
  real_T DerivativeGain_ha;            /* '<S211>/Derivative Gain' */
  real_T Filter_i;                     /* '<S211>/Filter' */
  real_T SumD_gv;                      /* '<S211>/SumD' */
  real_T FilterCoefficient_nn;         /* '<S211>/Filter Coefficient' */
  real_T Sum_bm;                       /* '<S211>/Sum' */
  real_T Sum23;                        /* '<S166>/Sum23' */
  real_T ProportionalGain_cc;          /* '<S210>/Proportional Gain' */
  real_T Integrator_fy;                /* '<S210>/Integrator' */
  real_T DerivativeGain_c;             /* '<S210>/Derivative Gain' */
  real_T Filter_g;                     /* '<S210>/Filter' */
  real_T SumD_b;                       /* '<S210>/SumD' */
  real_T FilterCoefficient_d;          /* '<S210>/Filter Coefficient' */
  real_T Sum_gh;                       /* '<S210>/Sum' */
  real_T ManualSwitch13;               /* '<S166>/Manual Switch13' */
  real_T deltafr_e;                    /* '<S166>/Gain6' */
  real_T deltafr_n;                    /* '<S166>/Manual Switch7' */
  real_T deltafl_e;                    /* '<S166>/Gain7' */
  real_T deltafl_f;                    /* '<S166>/Manual Switch8' */
  real_T ManualSwitch12;               /* '<S166>/Manual Switch12' */
  real_T ToBodyAxes4[3];               /* '<S166>/To Body Axes4' */
  real_T Sum14;                        /* '<S166>/Sum14' */
  real_T DerivativeGain_eb;            /* '<S200>/Derivative Gain' */
  real_T Filter_o0;                    /* '<S200>/Filter' */
  real_T SumD_o;                       /* '<S200>/SumD' */
  real_T FilterCoefficient_mw;         /* '<S200>/Filter Coefficient' */
  real_T IntegralGain_h;               /* '<S200>/Integral Gain' */
  real_T Integrator_m;                 /* '<S200>/Integrator' */
  real_T ProportionalGain_kp;          /* '<S200>/Proportional Gain' */
  real_T Sum_oz;                       /* '<S200>/Sum' */
  real_T delta;                        /* '<S166>/Sum2' */
  real_T Saturation9;                  /* '<S166>/Saturation9' */
  real_T ProportionalGain_n0;          /* '<S212>/Proportional Gain' */
  real_T Integrator_h;                 /* '<S212>/Integrator' */
  real_T DerivativeGain_b;             /* '<S212>/Derivative Gain' */
  real_T Filter_n;                     /* '<S212>/Filter' */
  real_T SumD_k;                       /* '<S212>/SumD' */
  real_T FilterCoefficient_p;          /* '<S212>/Filter Coefficient' */
  real_T Saturation3;                  /* '<S166>/Saturation3' */
  real_T Sum11;                        /* '<S166>/Sum11' */
  real_T ProportionalGain_i;           /* '<S202>/Proportional Gain' */
  real_T Integrator_fa;                /* '<S202>/Integrator' */
  real_T DerivativeGain_gz;            /* '<S202>/Derivative Gain' */
  real_T Filter_gk;                    /* '<S202>/Filter' */
  real_T SumD_kz;                      /* '<S202>/SumD' */
  real_T FilterCoefficient_b;          /* '<S202>/Filter Coefficient' */
  real_T ManualSwitch_o;               /* '<S166>/Manual Switch' */
  real_T ManualSwitch1_i;              /* '<S166>/Manual Switch1' */
  real_T ToBodyAxes6[3];               /* '<S166>/To Body Axes6' */
  real_T Sum8;                         /* '<S166>/Sum8' */
  real_T ProportionalGain_aq;          /* '<S213>/Proportional Gain' */
  real_T Integrator_j;                 /* '<S213>/Integrator' */
  real_T DerivativeGain_df;            /* '<S213>/Derivative Gain' */
  real_T Filter_if;                    /* '<S213>/Filter' */
  real_T SumD_ar;                      /* '<S213>/SumD' */
  real_T FilterCoefficient_kl;         /* '<S213>/Filter Coefficient' */
  real_T Sum_mf;                       /* '<S213>/Sum' */
  real_T Saturation2;                  /* '<S166>/Saturation2' */
  real_T ToBodyAxes7[3];               /* '<S166>/To Body Axes7' */
  real_T Sum3_b;                       /* '<S166>/Sum3' */
  real_T ProportionalGain_oc;          /* '<S205>/Proportional Gain' */
  real_T Integrator_lxr;               /* '<S205>/Integrator' */
  real_T DerivativeGain_m;             /* '<S205>/Derivative Gain' */
  real_T Filter_gu;                    /* '<S205>/Filter' */
  real_T SumD_hj;                      /* '<S205>/SumD' */
  real_T FilterCoefficient_i;          /* '<S205>/Filter Coefficient' */
  real_T Sum_lh;                       /* '<S205>/Sum' */
  real_T Saturation1;                  /* '<S166>/Saturation1' */
  real_T Sum6_n;                       /* '<S166>/Sum6' */
  real_T ProportionalGain_h;           /* '<S219>/Proportional Gain' */
  real_T Integrator_p;                 /* '<S219>/Integrator' */
  real_T DerivativeGain_k;             /* '<S219>/Derivative Gain' */
  real_T Filter_eg;                    /* '<S219>/Filter' */
  real_T SumD_fj;                      /* '<S219>/SumD' */
  real_T FilterCoefficient_mf;         /* '<S219>/Filter Coefficient' */
  real_T Sum_lq;                       /* '<S219>/Sum' */
  real_T Sum7_b;                       /* '<S166>/Sum7' */
  real_T ProportionalGain_le;          /* '<S215>/Proportional Gain' */
  real_T Integrator_gq;                /* '<S215>/Integrator' */
  real_T DerivativeGain_p;             /* '<S215>/Derivative Gain' */
  real_T Filter_ia;                    /* '<S215>/Filter' */
  real_T SumD_n;                       /* '<S215>/SumD' */
  real_T FilterCoefficient_ps;         /* '<S215>/Filter Coefficient' */
  real_T Sum_je;                       /* '<S215>/Sum' */
  real_T dRoll;                        /* '<S225>/dRoll' */
  real_T Gain_p;                       /* '<S225>/Gain' */
  real_T dThrottle;                    /* '<S223>/dThrottle' */
  real_T Bias_g;                       /* '<S223>/Bias' */
  real_T ManualSwitch6;                /* '<S166>/Manual Switch6' */
  real_T Sum9;                         /* '<S166>/Sum9' */
  real_T ProportionalGain_ky;          /* '<S220>/Proportional Gain' */
  real_T Integrator_a;                 /* '<S220>/Integrator' */
  real_T DerivativeGain_b1;            /* '<S220>/Derivative Gain' */
  real_T Filter_er;                    /* '<S220>/Filter' */
  real_T SumD_ae;                      /* '<S220>/SumD' */
  real_T FilterCoefficient_ip;         /* '<S220>/Filter Coefficient' */
  real_T Sum_asv;                      /* '<S220>/Sum' */
  real_T Sum10;                        /* '<S166>/Sum10' */
  real_T ProportionalGain_bf;          /* '<S216>/Proportional Gain' */
  real_T Integrator_bh;                /* '<S216>/Integrator' */
  real_T DerivativeGain_n5;            /* '<S216>/Derivative Gain' */
  real_T Filter_nz;                    /* '<S216>/Filter' */
  real_T SumD_bp;                      /* '<S216>/SumD' */
  real_T FilterCoefficient_pt;         /* '<S216>/Filter Coefficient' */
  real_T Sum_fu;                       /* '<S216>/Sum' */
  real_T dYaw;                         /* '<S226>/dYaw' */
  real_T Saturation7;                  /* '<S166>/Saturation7' */
  real_T ToBodyAxes5[3];               /* '<S166>/To Body Axes5' */
  real_T Sum15;                        /* '<S166>/Sum15' */
  real_T ProportionalGain_lw;          /* '<S204>/Proportional Gain' */
  real_T Integrator_av;                /* '<S204>/Integrator' */
  real_T DerivativeGain_o;             /* '<S204>/Derivative Gain' */
  real_T Filter_n4;                    /* '<S204>/Filter' */
  real_T SumD_f5;                      /* '<S204>/SumD' */
  real_T FilterCoefficient_o5;         /* '<S204>/Filter Coefficient' */
  real_T ManualSwitch9;                /* '<S166>/Manual Switch9' */
  real_T Sum12_m;                      /* '<S166>/Sum12' */
  real_T ProportionalGain_d;           /* '<S218>/Proportional Gain' */
  real_T Integrator_d2;                /* '<S218>/Integrator' */
  real_T DerivativeGain_au;            /* '<S218>/Derivative Gain' */
  real_T Filter_p;                     /* '<S218>/Filter' */
  real_T SumD_e;                       /* '<S218>/SumD' */
  real_T FilterCoefficient_m4;         /* '<S218>/Filter Coefficient' */
  real_T Sum_oy;                       /* '<S218>/Sum' */
  real_T PitchRateIn_g;                /* '<S166>/Sum13' */
  real_T ProportionalGain_at;          /* '<S214>/Proportional Gain' */
  real_T Integrator_bf;                /* '<S214>/Integrator' */
  real_T DerivativeGain_hn;            /* '<S214>/Derivative Gain' */
  real_T Filter_o0m;                   /* '<S214>/Filter' */
  real_T SumD_pe;                      /* '<S214>/SumD' */
  real_T FilterCoefficient_dq;         /* '<S214>/Filter Coefficient' */
  real_T Sum_mp;                       /* '<S214>/Sum' */
  real_T dPitch;                       /* '<S224>/dPitch' */
  real_T Throttle2noSaturation;        /* '<S221>/Sum2' */
  real_T Divide;                       /* '<S224>/Divide' */
  real_T Product_e4;                   /* '<S224>/Product' */
  real_T Gain_g;                       /* '<S226>/Gain' */
  real_T Gain1_n;                      /* '<S225>/Gain1' */
  real_T Throttle3noSaturation;        /* '<S221>/Sum1' */
  real_T Product1_b;                   /* '<S224>/Product1' */
  real_T Throttle4noSaturation;        /* '<S221>/Sum3' */
  real_T Gain1_b;                      /* '<S226>/Gain1' */
  real_T Throttle5noSaturation;        /* '<S221>/Sum' */
  real_T IntegralGain_m;               /* '<S201>/Integral Gain' */
  real_T IntegralGain_e;               /* '<S202>/Integral Gain' */
  real_T Sum_cs0;                      /* '<S166>/Sum' */
  real_T ProportionalGain_au;          /* '<S217>/Proportional Gain' */
  real_T Integrator_jk;                /* '<S217>/Integrator' */
  real_T DerivativeGain_l;             /* '<S217>/Derivative Gain' */
  real_T Filter_m;                     /* '<S217>/Filter' */
  real_T SumD_a1;                      /* '<S217>/SumD' */
  real_T FilterCoefficient_fx;         /* '<S217>/Filter Coefficient' */
  real_T Sum_i;                        /* '<S217>/Sum' */
  real_T Saturate_g;                   /* '<S217>/Saturate' */
  real_T Sum1_hp;                      /* '<S166>/Sum1' */
  real_T DerivativeGain_bm;            /* '<S203>/Derivative Gain' */
  real_T Filter_il;                    /* '<S203>/Filter' */
  real_T SumD_ng;                      /* '<S203>/SumD' */
  real_T FilterCoefficient_lp;         /* '<S203>/Filter Coefficient' */
  real_T IntegralGain_fs;              /* '<S203>/Integral Gain' */
  real_T IntegralGain_o;               /* '<S204>/Integral Gain' */
  real_T IntegralGain_g;               /* '<S205>/Integral Gain' */
  real_T IntegralGain_ay;              /* '<S206>/Integral Gain' */
  real_T IntegralGain_p;               /* '<S207>/Integral Gain' */
  real_T IntegralGain_k;               /* '<S208>/Integral Gain' */
  real_T IntegralGain_ox;              /* '<S209>/Integral Gain' */
  real_T IntegralGain_gz;              /* '<S210>/Integral Gain' */
  real_T IntegralGain_jv;              /* '<S211>/Integral Gain' */
  real_T IntegralGain_bz;              /* '<S212>/Integral Gain' */
  real_T IntegralGain_n;               /* '<S213>/Integral Gain' */
  real_T IntegralGain_d1;              /* '<S214>/Integral Gain' */
  real_T IntegralGain_ix;              /* '<S215>/Integral Gain' */
  real_T IntegralGain_kz;              /* '<S216>/Integral Gain' */
  real_T IntegralGain_b1;              /* '<S217>/Integral Gain' */
  real_T IntegralGain_da;              /* '<S218>/Integral Gain' */
  real_T IntegralGain_eu;              /* '<S219>/Integral Gain' */
  real_T IntegralGain_ao;              /* '<S220>/Integral Gain' */
  real_T Sum_mb;                       /* '<S212>/Sum' */
  real_T Sum_hg;                       /* '<S202>/Sum' */
  real_T Sum_ae;                       /* '<S208>/Sum' */
  real_T Sum25;                        /* '<S166>/Sum25' */
  real_T Gain8_l;                      /* '<S166>/Gain8' */
  real_T Sum17;                        /* '<S166>/Sum17' */
  real_T ManualSwitch11;               /* '<S166>/Manual Switch11' */
  real_T Sum_iw;                       /* '<S204>/Sum' */
  real_T Saturation6;                  /* '<S166>/Saturation6' */
  real_T TmpSignalConversionAtsincosInport1_i[3];
  real_T VectorConcatenate_lt[9];      /* '<S199>/Vector Concatenate' */
  real_T ToBodyAxes3_l[3];             /* '<S165>/To Body Axes3' */
  real_T ManualSwitch12_b;             /* '<S165>/Manual Switch12' */
  real_T ToBodyAxes4_p[3];             /* '<S165>/To Body Axes4' */
  real_T Sum14_o;                      /* '<S165>/Sum14' */
  real_T DerivativeGain_mg;            /* '<S179>/Derivative Gain' */
  real_T Filter_et;                    /* '<S179>/Filter' */
  real_T SumD_i0;                      /* '<S179>/SumD' */
  real_T FilterCoefficient_ct;         /* '<S179>/Filter Coefficient' */
  real_T IntegralGain_mk;              /* '<S179>/Integral Gain' */
  real_T Integrator_d5;                /* '<S179>/Integrator' */
  real_T ProportionalGain_hd;          /* '<S179>/Proportional Gain' */
  real_T Sum_pnb;                      /* '<S179>/Sum' */
  real_T delta_p;                      /* '<S165>/Sum2' */
  real_T Saturation9_k;                /* '<S165>/Saturation9' */
  real_T ProportionalGain_dp;          /* '<S184>/Proportional Gain' */
  real_T Integrator_n4;                /* '<S184>/Integrator' */
  real_T DerivativeGain_j4;            /* '<S184>/Derivative Gain' */
  real_T Filter_l3;                    /* '<S184>/Filter' */
  real_T SumD_j;                       /* '<S184>/SumD' */
  real_T FilterCoefficient_dn;         /* '<S184>/Filter Coefficient' */
  real_T Saturation3_o;                /* '<S165>/Saturation3' */
  real_T Sum11_j;                      /* '<S165>/Sum11' */
  real_T ProportionalGain_lh;          /* '<S180>/Proportional Gain' */
  real_T Integrator_bx;                /* '<S180>/Integrator' */
  real_T DerivativeGain_in;            /* '<S180>/Derivative Gain' */
  real_T Filter_fe;                    /* '<S180>/Filter' */
  real_T SumD_ji;                      /* '<S180>/SumD' */
  real_T FilterCoefficient_im;         /* '<S180>/Filter Coefficient' */
  real_T ManualSwitch_k;               /* '<S165>/Manual Switch' */
  real_T ManualSwitch1_o;              /* '<S165>/Manual Switch1' */
  real_T ToBodyAxes6_a[3];             /* '<S165>/To Body Axes6' */
  real_T Sum8_h;                       /* '<S165>/Sum8' */
  real_T ProportionalGain_m;           /* '<S185>/Proportional Gain' */
  real_T Integrator_iz;                /* '<S185>/Integrator' */
  real_T DerivativeGain_j3;            /* '<S185>/Derivative Gain' */
  real_T Filter_j;                     /* '<S185>/Filter' */
  real_T SumD_g3;                      /* '<S185>/SumD' */
  real_T FilterCoefficient_dj;         /* '<S185>/Filter Coefficient' */
  real_T Sum_pnh;                      /* '<S185>/Sum' */
  real_T Saturation2_p;                /* '<S165>/Saturation2' */
  real_T ToBodyAxes7_e[3];             /* '<S165>/To Body Axes7' */
  real_T Sum3_g;                       /* '<S165>/Sum3' */
  real_T ProportionalGain_m0;          /* '<S183>/Proportional Gain' */
  real_T Integrator_fay;               /* '<S183>/Integrator' */
  real_T DerivativeGain_je;            /* '<S183>/Derivative Gain' */
  real_T Filter_a;                     /* '<S183>/Filter' */
  real_T SumD_nt;                      /* '<S183>/SumD' */
  real_T FilterCoefficient_dr;         /* '<S183>/Filter Coefficient' */
  real_T Sum_fr;                       /* '<S183>/Sum' */
  real_T Saturation1_h;                /* '<S165>/Saturation1' */
  real_T Sum6_p;                       /* '<S165>/Sum6' */
  real_T ProportionalGain_iq;          /* '<S191>/Proportional Gain' */
  real_T Integrator_k;                 /* '<S191>/Integrator' */
  real_T DerivativeGain_a1;            /* '<S191>/Derivative Gain' */
  real_T Filter_fd;                    /* '<S191>/Filter' */
  real_T SumD_e0;                      /* '<S191>/SumD' */
  real_T FilterCoefficient_h;          /* '<S191>/Filter Coefficient' */
  real_T Sum_a2b;                      /* '<S191>/Sum' */
  real_T Sum7_l;                       /* '<S165>/Sum7' */
  real_T ProportionalGain_dv;          /* '<S187>/Proportional Gain' */
  real_T Integrator_af;                /* '<S187>/Integrator' */
  real_T DerivativeGain_na;            /* '<S187>/Derivative Gain' */
  real_T Filter_k;                     /* '<S187>/Filter' */
  real_T SumD_m;                       /* '<S187>/SumD' */
  real_T FilterCoefficient_dv;         /* '<S187>/Filter Coefficient' */
  real_T Sum_af;                       /* '<S187>/Sum' */
  real_T dRoll_m;                      /* '<S197>/dRoll' */
  real_T Gain_h;                       /* '<S197>/Gain' */
  real_T dThrottle_g;                  /* '<S195>/dThrottle' */
  real_T Bias_f;                       /* '<S195>/Bias' */
  real_T ManualSwitch6_e;              /* '<S165>/Manual Switch6' */
  real_T Sum9_e;                       /* '<S165>/Sum9' */
  real_T ProportionalGain_nl;          /* '<S192>/Proportional Gain' */
  real_T Integrator_ar;                /* '<S192>/Integrator' */
  real_T DerivativeGain_oq;            /* '<S192>/Derivative Gain' */
  real_T Filter_b;                     /* '<S192>/Filter' */
  real_T SumD_c2;                      /* '<S192>/SumD' */
  real_T FilterCoefficient_cz;         /* '<S192>/Filter Coefficient' */
  real_T Sum_i4;                       /* '<S192>/Sum' */
  real_T Sum10_n;                      /* '<S165>/Sum10' */
  real_T ProportionalGain_p;           /* '<S188>/Proportional Gain' */
  real_T Integrator_ly;                /* '<S188>/Integrator' */
  real_T DerivativeGain_li;            /* '<S188>/Derivative Gain' */
  real_T Filter_p0;                    /* '<S188>/Filter' */
  real_T SumD_lv;                      /* '<S188>/SumD' */
  real_T FilterCoefficient_ix;         /* '<S188>/Filter Coefficient' */
  real_T Sum_dv;                       /* '<S188>/Sum' */
  real_T dYaw_a;                       /* '<S198>/dYaw' */
  real_T Saturation7_k;                /* '<S165>/Saturation7' */
  real_T ToBodyAxes5_b[3];             /* '<S165>/To Body Axes5' */
  real_T Sum15_g;                      /* '<S165>/Sum15' */
  real_T ProportionalGain_fn;          /* '<S182>/Proportional Gain' */
  real_T Integrator_ne;                /* '<S182>/Integrator' */
  real_T DerivativeGain_fe;            /* '<S182>/Derivative Gain' */
  real_T Filter_ik;                    /* '<S182>/Filter' */
  real_T SumD_mm;                      /* '<S182>/SumD' */
  real_T FilterCoefficient_g;          /* '<S182>/Filter Coefficient' */
  real_T Sum_ib;                       /* '<S182>/Sum' */
  real_T Saturation6_k;                /* '<S165>/Saturation6' */
  real_T Sum12_e;                      /* '<S165>/Sum12' */
  real_T ProportionalGain_dx;          /* '<S190>/Proportional Gain' */
  real_T Integrator_nen;               /* '<S190>/Integrator' */
  real_T DerivativeGain_no;            /* '<S190>/Derivative Gain' */
  real_T Filter_h5;                    /* '<S190>/Filter' */
  real_T SumD_fi;                      /* '<S190>/SumD' */
  real_T FilterCoefficient_dc;         /* '<S190>/Filter Coefficient' */
  real_T Sum_cy;                       /* '<S190>/Sum' */
  real_T PitchRateIn_c;                /* '<S165>/Sum13' */
  real_T ProportionalGain_ab;          /* '<S186>/Proportional Gain' */
  real_T Integrator_hd;                /* '<S186>/Integrator' */
  real_T DerivativeGain_cp;            /* '<S186>/Derivative Gain' */
  real_T Filter_dy;                    /* '<S186>/Filter' */
  real_T SumD_ei;                      /* '<S186>/SumD' */
  real_T FilterCoefficient_mr;         /* '<S186>/Filter Coefficient' */
  real_T Sum_pf;                       /* '<S186>/Sum' */
  real_T dPitch_h;                     /* '<S196>/dPitch' */
  real_T Throttle2noSaturation_l;      /* '<S193>/Sum2' */
  real_T Divide_h;                     /* '<S196>/Divide' */
  real_T Product_ma;                   /* '<S196>/Product' */
  real_T Gain_ac;                      /* '<S198>/Gain' */
  real_T Gain1_h;                      /* '<S197>/Gain1' */
  real_T Throttle3noSaturation_d;      /* '<S193>/Sum1' */
  real_T Product1_l;                   /* '<S196>/Product1' */
  real_T Throttle4noSaturation_d;      /* '<S193>/Sum3' */
  real_T Gain1_ne;                     /* '<S198>/Gain1' */
  real_T Throttle5noSaturation_g;      /* '<S193>/Sum' */
  real_T IntegralGain_nx;              /* '<S180>/Integral Gain' */
  real_T Sum_ly;                       /* '<S165>/Sum' */
  real_T ProportionalGain_nb;          /* '<S189>/Proportional Gain' */
  real_T Integrator_gb;                /* '<S189>/Integrator' */
  real_T DerivativeGain_j33;           /* '<S189>/Derivative Gain' */
  real_T Filter_kh;                    /* '<S189>/Filter' */
  real_T SumD_c3;                      /* '<S189>/SumD' */
  real_T FilterCoefficient_ko;         /* '<S189>/Filter Coefficient' */
  real_T Sum_kb;                       /* '<S189>/Sum' */
  real_T Saturate_iq;                  /* '<S189>/Saturate' */
  real_T Sum1_c;                       /* '<S165>/Sum1' */
  real_T DerivativeGain_inf;           /* '<S181>/Derivative Gain' */
  real_T Filter_gp;                    /* '<S181>/Filter' */
  real_T SumD_l3;                      /* '<S181>/SumD' */
  real_T FilterCoefficient_lz;         /* '<S181>/Filter Coefficient' */
  real_T IntegralGain_ke;              /* '<S181>/Integral Gain' */
  real_T IntegralGain_l;               /* '<S182>/Integral Gain' */
  real_T IntegralGain_g4;              /* '<S183>/Integral Gain' */
  real_T IntegralGain_hx;              /* '<S184>/Integral Gain' */
  real_T IntegralGain_c;               /* '<S185>/Integral Gain' */
  real_T IntegralGain_ho;              /* '<S186>/Integral Gain' */
  real_T IntegralGain_kc;              /* '<S187>/Integral Gain' */
  real_T IntegralGain_m1;              /* '<S188>/Integral Gain' */
  real_T IntegralGain_df;              /* '<S189>/Integral Gain' */
  real_T IntegralGain_an;              /* '<S190>/Integral Gain' */
  real_T IntegralGain_pv;              /* '<S191>/Integral Gain' */
  real_T IntegralGain_pn;              /* '<S192>/Integral Gain' */
  real_T Sum_dk;                       /* '<S184>/Sum' */
  real_T Sum_ne;                       /* '<S180>/Sum' */
  real_T ToBodyAxes1[3];               /* '<S164>/To Body Axes1' */
  real_T VbodyIn;                      /* '<S164>/Sum7' */
  real_T ProportionalGain_a3;          /* '<S176>/Proportional Gain' */
  real_T Integrator_ms;                /* '<S176>/Integrator' */
  real_T DerivativeGain_nt;            /* '<S176>/Derivative Gain' */
  real_T Filter_mm;                    /* '<S176>/Filter' */
  real_T SumD_gx;                      /* '<S176>/SumD' */
  real_T FilterCoefficient_ge;         /* '<S176>/Filter Coefficient' */
  real_T Sum_bi;                       /* '<S176>/Sum' */
  real_T Saturate_j;                   /* '<S176>/Saturate' */
  real_T Sum1_k;                       /* '<S164>/Sum1' */
  real_T Saturation_m;                 /* '<S164>/Saturation' */
  real_T ProportionalGain_c1;          /* '<S169>/Proportional Gain' */
  real_T Integrator_la;                /* '<S169>/Integrator' */
  real_T DerivativeGain_kp;            /* '<S169>/Derivative Gain' */
  real_T Filter_m2;                    /* '<S169>/Filter' */
  real_T SumD_me;                      /* '<S169>/SumD' */
  real_T FilterCoefficient_f0;         /* '<S169>/Filter Coefficient' */
  real_T Sum_gg;                       /* '<S169>/Sum' */
  real_T Saturate_o;                   /* '<S169>/Saturate' */
  real_T ManualSwitch_d;               /* '<S164>/Manual Switch' */
  real_T Sum12_o;                      /* '<S164>/Sum12' */
  real_T ProportionalGain_cd;          /* '<S171>/Proportional Gain' */
  real_T Integrator_ew;                /* '<S171>/Integrator' */
  real_T DerivativeGain_fn;            /* '<S171>/Derivative Gain' */
  real_T Filter_b4;                    /* '<S171>/Filter' */
  real_T SumD_fjm;                     /* '<S171>/SumD' */
  real_T FilterCoefficient_n0;         /* '<S171>/Filter Coefficient' */
  real_T Sum_j5;                       /* '<S171>/Sum' */
  real_T PitchRateIn_h;                /* '<S164>/Sum13' */
  real_T ProportionalGain_j;           /* '<S170>/Proportional Gain' */
  real_T Integrator_lp4;               /* '<S170>/Integrator' */
  real_T DerivativeGain_m4;            /* '<S170>/Derivative Gain' */
  real_T Filter_fn;                    /* '<S170>/Filter' */
  real_T SumD_na;                      /* '<S170>/SumD' */
  real_T FilterCoefficient_jh;         /* '<S170>/Filter Coefficient' */
  real_T Sum_it;                       /* '<S170>/Sum' */
  real_T Saturate_ia;                  /* '<S170>/Saturate' */
  real_T deltae_f;                     /* '<S164>/Gain4' */
  real_T Sum3_i;                       /* '<S164>/Sum3' */
  real_T ProportionalGain_i2;          /* '<S173>/Proportional Gain' */
  real_T Integrator_et;                /* '<S173>/Integrator' */
  real_T DerivativeGain_me;            /* '<S173>/Derivative Gain' */
  real_T Filter_p4;                    /* '<S173>/Filter' */
  real_T SumD_hb;                      /* '<S173>/SumD' */
  real_T FilterCoefficient_d3;         /* '<S173>/Filter Coefficient' */
  real_T Sum_i1;                       /* '<S173>/Sum' */
  real_T Sum2_aq;                      /* '<S164>/Sum2' */
  real_T ProportionalGain_ce;          /* '<S172>/Proportional Gain' */
  real_T Integrator_o;                 /* '<S172>/Integrator' */
  real_T DerivativeGain_jw;            /* '<S172>/Derivative Gain' */
  real_T Filter_c;                     /* '<S172>/Filter' */
  real_T SumD_jm;                      /* '<S172>/SumD' */
  real_T FilterCoefficient_hr;         /* '<S172>/Filter Coefficient' */
  real_T deltar_k;                     /* '<S164>/Gain5' */
  real_T Sum5_c;                       /* '<S164>/Sum5' */
  real_T ProportionalGain_bl;          /* '<S175>/Proportional Gain' */
  real_T Integrator_lr;                /* '<S175>/Integrator' */
  real_T DerivativeGain_hj;            /* '<S175>/Derivative Gain' */
  real_T Filter_dp;                    /* '<S175>/Filter' */
  real_T SumD_kp;                      /* '<S175>/SumD' */
  real_T FilterCoefficient_ok;         /* '<S175>/Filter Coefficient' */
  real_T Sum_mm;                       /* '<S175>/Sum' */
  real_T Sum4_k;                       /* '<S164>/Sum4' */
  real_T ProportionalGain_j2;          /* '<S174>/Proportional Gain' */
  real_T Integrator_ey;                /* '<S174>/Integrator' */
  real_T DerivativeGain_ic;            /* '<S174>/Derivative Gain' */
  real_T Filter_eo;                    /* '<S174>/Filter' */
  real_T SumD_pf;                      /* '<S174>/SumD' */
  real_T FilterCoefficient_cv;         /* '<S174>/Filter Coefficient' */
  real_T Sum_gf;                       /* '<S174>/Sum' */
  real_T deltafr_c;                    /* '<S164>/Gain6' */
  real_T deltafl_p;                    /* '<S164>/Gain7' */
  real_T IntegralGain_m3;              /* '<S169>/Integral Gain' */
  real_T IntegralGain_jo;              /* '<S170>/Integral Gain' */
  real_T IntegralGain_ah;              /* '<S171>/Integral Gain' */
  real_T IntegralGain_fk;              /* '<S172>/Integral Gain' */
  real_T IntegralGain_pd;              /* '<S173>/Integral Gain' */
  real_T IntegralGain_hb;              /* '<S174>/Integral Gain' */
  real_T IntegralGain_a3;              /* '<S175>/Integral Gain' */
  real_T IntegralGain_cv;              /* '<S176>/Integral Gain' */
  real_T ManualSwitch7;                /* '<S164>/Manual Switch7' */
  real_T Sum_ggh;                      /* '<S172>/Sum' */
  real_T Sum6_g;                       /* '<S164>/Sum6' */
  real_T Gain8_i;                      /* '<S164>/Gain8' */
  real_T Sum_ast;                      /* '<S164>/Sum' */
  real_T Product_jq;                   /* '<S158>/Product' */
  real_T sp2;                          /* '<S108>/sp[2]' */
  real_T cp2;                          /* '<S108>/cp[2]' */
  real_T Gain_k[11];                   /* '<S108>/Gain' */
  real_T Gain1_e[11];                  /* '<S108>/Gain1' */
  real_T OutportBufferForcp13[13];
  real_T OutportBufferForsp13[13];
  real_T b2;                           /* '<S107>/b2' */
  real_T a2;                           /* '<S107>/a2' */
  real_T c2;                           /* '<S150>/Sum1' */
  real_T Product_gg;                   /* '<S150>/Product' */
  real_T Sum_l4;                       /* '<S150>/Sum' */
  real_T q1;                           /* '<S107>/Product1' */
  real_T Product9_i;                   /* '<S149>/Product9' */
  real_T Product10_o;                  /* '<S149>/Product10' */
  real_T Sum7_g;                       /* '<S149>/Sum7' */
  real_T Sum8_p;                       /* '<S147>/Sum8' */
  real_T Gain_f;                       /* '<S152>/Gain' */
  real_T Product6_p;                   /* '<S152>/Product6' */
  real_T a4;                           /* '<S152>/a4' */
  real_T c4;                           /* '<S152>/Sum9' */
  real_T Product8_d;                   /* '<S152>/Product8' */
  real_T Sum6_j;                       /* '<S152>/Sum6' */
  real_T Product1_oa;                  /* '<S152>/Product1' */
  real_T Product7_a;                   /* '<S152>/Product7' */
  real_T r2;                           /* '<S152>/Sum5' */
  real_T sqrt_l;                       /* '<S107>/sqrt' */
  real_T Product11_j;                  /* '<S147>/Product11' */
  real_T Sum2_if;                      /* '<S151>/Sum2' */
  real_T Sum1_n;                       /* '<S151>/Sum1' */
  real_T Product1_j;                   /* '<S151>/Product1' */
  real_T q2;                           /* '<S151>/Product2' */
  real_T Product3_c;                   /* '<S148>/Product3' */
  real_T Sum3_c;                       /* '<S148>/Sum3' */
  real_T Product4_l;                   /* '<S148>/Product4' */
  real_T Product1_m4;                  /* '<S153>/Product1' */
  real_T c2_e;                         /* '<S153>/Sum1' */
  real_T Product12;                    /* '<S153>/Product12' */
  real_T Product5_h;                   /* '<S154>/Product5' */
  real_T Sum4_g;                       /* '<S154>/Sum4' */
  real_T sqrt_m;                       /* '<S154>/sqrt' */
  real_T Sum1_lw[4];                   /* '<S106>/Sum1' */
  real_T Sum1_bu;                      /* '<S115>/Sum1' */
  real_T Sum2_p;                       /* '<S115>/Sum2' */
  real_T Sum3_co;                      /* '<S115>/Sum3' */
  real_T Sum5_k;                       /* '<S115>/Sum5' */
  real_T Merge_l[169];                 /* '<S144>/Merge' */
  real_T Sum2_d[169];                  /* '<S117>/Sum2' */
  real_T Merge1_j;                     /* '<S116>/Merge1' */
  real_T Assignment_i[169];            /* '<S116>/Assignment' */
  real_T Merge_im;                     /* '<S116>/Merge' */
  real_T Assignment_snorm[169];        /* '<S116>/Assignment_snorm' */
  real_T TmpSignalConversionAtppnInport1[13];
  real_T Product2_o;                   /* '<S118>/Product2' */
  real_T Assignment2_o[13];            /* '<S123>/Assignment2' */
  real_T Assignment2_e[13];            /* '<S122>/Assignment2' */
  real_T Product2_i[2];                /* '<S88>/Product2' */
  real_T Product1_mk[2];               /* '<S88>/Product1' */
  real_T VectorConcatenate_j3[3];      /* '<S87>/Vector Concatenate' */
  real_T Product_lt[3];                /* '<S87>/Product' */
  real_T Sum2_g[3];                    /* '<S83>/Sum2' */
  real_T Sum1_jc;                      /* '<S83>/Sum1' */
  real_T Sum_em;                       /* '<S83>/Sum' */
  real_T Product1_hj[3];               /* '<S83>/Product1' */
  real_T Product2_oo[2];               /* '<S90>/Product2' */
  real_T Product1_o0[2];               /* '<S90>/Product1' */
  real_T VectorConcatenate_a[3];       /* '<S89>/Vector Concatenate' */
  real_T Product_mb[3];                /* '<S89>/Product' */
  real_T Product2_my[2];               /* '<S80>/Product2' */
  real_T Product1_k0[2];               /* '<S80>/Product1' */
  real_T VectorConcatenate_p[3];       /* '<S79>/Vector Concatenate' */
  real_T Product_gv[3];                /* '<S79>/Product' */
  real_T Sum2_g5[3];                   /* '<S75>/Sum2' */
  real_T Sum1_jk;                      /* '<S75>/Sum1' */
  real_T Sum_eg;                       /* '<S75>/Sum' */
  real_T Product1_h5[3];               /* '<S75>/Product1' */
  real_T Product2_id[2];               /* '<S82>/Product2' */
  real_T Product1_dl[2];               /* '<S82>/Product1' */
  real_T VectorConcatenate_aw[3];      /* '<S81>/Vector Concatenate' */
  real_T Product_kg[3];                /* '<S81>/Product' */
  real_T VLwg[2];                      /* '<S72>/V//Lwg' */
  real_T u[2];                         /* '<S72>/2' */
  real_T LugV1[2];                     /* '<S72>/Lug//V1' */
  real_T dt[2];                        /* '<S72>/dt' */
  real_T Sum1_a[2];                    /* '<S72>/Sum1' */
  real_T UnitDelay[2];                 /* '<S72>/Unit Delay' */
  real_T LugV2[2];                     /* '<S72>/Lug//V2' */
  real_T Sum_da[2];                    /* '<S72>/Sum' */
  real_T VLvg[2];                      /* '<S71>/V//Lvg' */
  real_T u_h[2];                       /* '<S71>/2' */
  real_T LugV1_m[2];                   /* '<S71>/Lug//V1' */
  real_T dt_g[2];                      /* '<S71>/dt' */
  real_T Sum1_bg[2];                   /* '<S71>/Sum1' */
  real_T UnitDelay_d[2];               /* '<S71>/Unit Delay' */
  real_T LugV2_d[2];                   /* '<S71>/Lug//V2' */
  real_T Sum_nt[2];                    /* '<S71>/Sum' */
  real_T VLug[2];                      /* '<S70>/V//Lug' */
  real_T u_i[2];                       /* '<S70>/2' */
  real_T LugV1_l[2];                   /* '<S70>/Lug//V1' */
  real_T dt_l[2];                      /* '<S70>/dt' */
  real_T Sum1_o5[2];                   /* '<S70>/Sum1' */
  real_T UnitDelay_k[2];               /* '<S70>/Unit Delay' */
  real_T LugV2_c[2];                   /* '<S70>/Lug//V2' */
  real_T Sum_eq[2];                    /* '<S70>/Sum' */
  real_T dt1;                          /* '<S69>/dt1' */
  real_T ar;                           /* '<S69>/w1' */
  real_T dt_i;                         /* '<S69>/dt' */
  real_T Sum2_fh;                      /* '<S69>/Sum2' */
  real_T UnitDelay_b[2];               /* '<S69>/Unit Delay' */
  real_T LugV2_p[2];                   /* '<S69>/Lug//V2' */
  real_T UnitDelay1[2];                /* '<S69>/Unit Delay1' */
  real_T Sum3_jj[2];                   /* '<S69>/Sum3' */
  real_T w2[2];                        /* '<S69>/w2' */
  real_T Sum1_p2[2];                   /* '<S69>/Sum1' */
  real_T dt1_k;                        /* '<S68>/dt1' */
  real_T aq;                           /* '<S68>/w1' */
  real_T dt_f;                         /* '<S68>/dt' */
  real_T Sum2_l;                       /* '<S68>/Sum2' */
  real_T UnitDelay_n[2];               /* '<S68>/Unit Delay' */
  real_T LugV2_c3[2];                  /* '<S68>/Lug//V2' */
  real_T UnitDelay1_j[2];              /* '<S68>/Unit Delay1' */
  real_T Sum3_p[2];                    /* '<S68>/Sum3' */
  real_T w2_b[2];                      /* '<S68>/w2' */
  real_T Sum1_gx[2];                   /* '<S68>/Sum1' */
  real_T w2_bi[2];                     /* '<S67>/w2' */
  real_T ap[2];                        /* '<S67>/w1' */
  real_T u_j[2];                       /* '<S67>/2' */
  real_T w4[2];                        /* '<S67>/w4' */
  real_T sp[2];                        /* '<S67>/w3' */
  real_T LugV1_f[2];                   /* '<S67>/Lug//V1' */
  real_T dt_d[2];                      /* '<S67>/dt' */
  real_T Sum1_cw[2];                   /* '<S67>/Sum1' */
  real_T UnitDelay_o[2];               /* '<S67>/Unit Delay' */
  real_T LugV2_pf[2];                  /* '<S67>/Lug//V2' */
  real_T Sum_ax[2];                    /* '<S67>/Sum' */
  real_T DistanceintoGustxLimitedtogustlengthd;/* '<S50>/Distance into Gust (x) (Limited to gust length d)' */
  real_T uftinf;                       /* '<S49>/3ft-->inf' */
  real_T hz0;                          /* '<S49>/h//z0' */
  real_T Product_dl;                   /* '<S49>/Product' */
  real_T Product1_ek[3];               /* '<S49>/Product1' */
  real_T TransformfromInertialtoBodyaxes[3];/* '<S49>/Transform from Inertial to Body axes' */
  real_T WindVelocity[3];              /* '<S15>/Sum' */
  real_T Throttle1;                    /* '<S163>/IC' */
  real_T Throttle2_g;                  /* '<S163>/IC1' */
  real_T Throttle3_jz;                 /* '<S163>/IC2' */
  real_T Throttle4_j;                  /* '<S163>/IC3' */
  real_T Throttle5_a;                  /* '<S163>/IC4' */
  int32_T DataTypeConversion;          /* '<S167>/Data Type Conversion' */
  uint32_T PreLookUpIndexSearchaltitude_o1;/* '<S73>/PreLook-Up Index Search  (altitude)' */
  uint32_T PreLookUpIndexSearchprobofexceed_o1;/* '<S73>/PreLook-Up Index Search  (prob of exceed)' */
  uint32_T VoltagePrelookup_o1;        /* '<S353>/Voltage Prelookup' */
  uint32_T HeightPrelookup_o1;         /* '<S353>/Height Prelookup' */
  uint32_T VflightPrelookup_o1;        /* '<S353>/Vflight Prelookup' */
  uint32_T VoltagePrelookup_o1_a;      /* '<S367>/Voltage Prelookup' */
  uint32_T HeightPrelookup_o1_g;       /* '<S367>/Height Prelookup' */
  uint32_T VflightPrelookup_o1_e;      /* '<S367>/Vflight Prelookup' */
  uint32_T VoltagePrelookup_o1_l;      /* '<S381>/Voltage Prelookup' */
  uint32_T HeightPrelookup_o1_i;       /* '<S381>/Height Prelookup' */
  uint32_T VflightPrelookup_o1_m;      /* '<S381>/Vflight Prelookup' */
  uint32_T VoltagePrelookup_o1_i;      /* '<S395>/Voltage Prelookup' */
  uint32_T HeightPrelookup_o1_d;       /* '<S395>/Height Prelookup' */
  uint32_T VflightPrelookup_o1_a;      /* '<S395>/Vflight Prelookup' */
  uint32_T VoltagePrelookup_o1_a2;     /* '<S409>/Voltage Prelookup' */
  uint32_T HeightPrelookup_o1_l;       /* '<S409>/Height Prelookup' */
  uint32_T VflightPrelookup_o1_g;      /* '<S409>/Vflight Prelookup' */
  uint32_T VflightPrelookup_o1_d;      /* '<S410>/Vflight Prelookup' */
  uint32_T RPMPrelookup_o1;            /* '<S410>/RPM Prelookup' */
  uint32_T VflightPrelookup_o1_h;      /* '<S396>/Vflight Prelookup' */
  uint32_T RPMPrelookup_o1_k;          /* '<S396>/RPM Prelookup' */
  uint32_T VflightPrelookup_o1_gx;     /* '<S382>/Vflight Prelookup' */
  uint32_T RPMPrelookup_o1_c;          /* '<S382>/RPM Prelookup' */
  uint32_T VflightPrelookup_o1_j;      /* '<S368>/Vflight Prelookup' */
  uint32_T RPMPrelookup_o1_cy;         /* '<S368>/RPM Prelookup' */
  uint32_T VflightPrelookup_o1_f;      /* '<S354>/Vflight Prelookup' */
  uint32_T RPMPrelookup_o1_f;          /* '<S354>/RPM Prelookup' */
  boolean_T RelationalOperator;        /* '<S112>/Relational Operator' */
  boolean_T RelationalOperator_b;      /* '<S111>/Relational Operator' */
  boolean_T
    HiddenBuf_InsertedFor_Convertfromgeodetictosphericalcoordinates_at_inport_2;/* '<S104>/Has longitude changed ' */
  boolean_T LogicalOperator;           /* '<S110>/Logical Operator' */
  boolean_T
    HiddenBuf_InsertedFor_Convertfromgeodetictosphericalcoordinates_at_inport_5;/* '<S104>/Has altitude or latitude changed' */
  boolean_T LogicalOperator2;          /* '<S47>/Logical Operator2' */
  boolean_T HiddenBuf_InsertedFor_Distanceintogustx_at_inport_1;/* '<S47>/Logical Operator2' */
  boolean_T LogicalOperator1;          /* '<S47>/Logical Operator1' */
  boolean_T HiddenBuf_InsertedFor_Distanceintogusty_at_inport_1;/* '<S47>/Logical Operator1' */
  boolean_T LogicalOperator3;          /* '<S47>/Logical Operator3' */
  boolean_T HiddenBuf_InsertedFor_Distanceintogustz_at_inport_1;/* '<S47>/Logical Operator3' */
  boolean_T conjunction;               /* '<S96>/conjunction' */
  boolean_T conjunction_g;             /* '<S97>/conjunction' */
  boolean_T conjunction_d;             /* '<S98>/conjunction' */
  boolean_T conjunction_k;             /* '<S100>/conjunction' */
  boolean_T buttons[12];               /* '<S5>/Joystick Input' */
  boolean_T Compare;                   /* '<S228>/Compare' */
  boolean_T HiddenBuf_InsertedFor_Quadcopter_at_inport_2;/* '<S3>/Select Active FCS' */
  boolean_T Compare_j;                 /* '<S229>/Compare' */
  boolean_T HiddenBuf_InsertedFor_QuadcopterFixedWing_at_inport_2;/* '<S3>/Select Active FCS' */
  boolean_T Compare_g;                 /* '<S230>/Compare' */
  boolean_T HiddenBuf_InsertedFor_FixedWingClimb_at_inport_2;/* '<S3>/Select Active FCS' */
  boolean_T Compare_a;                 /* '<S232>/Compare' */
  boolean_T HiddenBuf_InsertedFor_TakeOff_at_inport_2;/* '<S3>/Select Active FCS' */
  B_Distanceintogusty_mainV03_56_T Distanceintogustz;/* '<S47>/Distance into gust (z)' */
  B_Distanceintogusty_mainV03_56_T Distanceintogusty;/* '<S47>/Distance into gust (y)' */
} B_mainV03_56_T;

/* Block states (auto storage) for system '<Root>' */
typedef struct {
  real_T UnitDelay_DSTATE;             /* '<S464>/Unit Delay' */
  real_T UnitDelay1_DSTATE[2];         /* '<S155>/Unit Delay1' */
  real_T UnitDelay1_DSTATE_m;          /* '<S115>/Unit Delay1' */
  real_T UnitDelay3_DSTATE;            /* '<S115>/Unit Delay3' */
  real_T UnitDelay2_DSTATE_i;          /* '<S115>/Unit Delay2' */
  real_T UnitDelay4_DSTATE;            /* '<S115>/Unit Delay4' */
  real_T UnitDelay_DSTATE_f[169];      /* '<S117>/Unit Delay' */
  real_T UnitDelay_DSTATE_i[169];      /* '<S144>/Unit Delay' */
  real_T UnitDelay_DSTATE_k[169];      /* '<S116>/Unit Delay' */
  real_T UnitDelay1_DSTATE_k[169];     /* '<S116>/Unit Delay1' */
  real_T UnitDelay1_DSTATE_n[13];      /* '<S118>/Unit Delay1' */
  real_T UnitDelay_DSTATE_l[2];        /* '<S72>/Unit Delay' */
  real_T UnitDelay_DSTATE_o[2];        /* '<S71>/Unit Delay' */
  real_T UnitDelay_DSTATE_p[2];        /* '<S70>/Unit Delay' */
  real_T UnitDelay_DSTATE_ju[2];       /* '<S69>/Unit Delay' */
  real_T UnitDelay1_DSTATE_ml[2];      /* '<S69>/Unit Delay1' */
  real_T UnitDelay_DSTATE_fj[2];       /* '<S68>/Unit Delay' */
  real_T UnitDelay1_DSTATE_j[2];       /* '<S68>/Unit Delay1' */
  real_T UnitDelay_DSTATE_c[2];        /* '<S67>/Unit Delay' */
  real_T IC4_FirstOutputTime;          /* '<S4>/IC4' */
  real_T IC11_FirstOutputTime;         /* '<S4>/IC11' */
  real_T IC12_FirstOutputTime;         /* '<S4>/IC12' */
  real_T IC2_FirstOutputTime;          /* '<S4>/IC2' */
  real_T IC5_FirstOutputTime;          /* '<S4>/IC5' */
  real_T IC3_FirstOutputTime;          /* '<S4>/IC3' */
  real_T IC7_FirstOutputTime;          /* '<S4>/IC7' */
  real_T IC8_FirstOutputTime;          /* '<S4>/IC8' */
  real_T TimeStampA;                   /* '<S4>/Derivative' */
  real_T LastUAtTimeA;                 /* '<S4>/Derivative' */
  real_T TimeStampB;                   /* '<S4>/Derivative' */
  real_T LastUAtTimeB;                 /* '<S4>/Derivative' */
  real_T PrevYA;                       /* '<S4>/Rate Limiter' */
  real_T PrevYB;                       /* '<S4>/Rate Limiter' */
  real_T LastMajorTimeA;               /* '<S4>/Rate Limiter' */
  real_T LastMajorTimeB;               /* '<S4>/Rate Limiter' */
  real_T TimeStampA_i;                 /* '<S4>/Derivative1' */
  real_T LastUAtTimeA_o;               /* '<S4>/Derivative1' */
  real_T TimeStampB_a;                 /* '<S4>/Derivative1' */
  real_T LastUAtTimeB_j;               /* '<S4>/Derivative1' */
  real_T PrevYA_o;                     /* '<S4>/Rate Limiter1' */
  real_T PrevYB_d;                     /* '<S4>/Rate Limiter1' */
  real_T LastMajorTimeA_f;             /* '<S4>/Rate Limiter1' */
  real_T LastMajorTimeB_d;             /* '<S4>/Rate Limiter1' */
  real_T SFunction_temp_table[8];      /* '<S12>/S-Function' */
  real_T SFunction_pres_table[8];      /* '<S12>/S-Function' */
  real_T WGS84GravitySFunction_h;      /* '<S14>/WGS84 Gravity S-Function' */
  real_T WGS84GravitySFunction_phi;    /* '<S14>/WGS84 Gravity S-Function' */
  real_T WGS84GravitySFunction_lambda; /* '<S14>/WGS84 Gravity S-Function' */
  real_T WGS84GravitySFunction_gamma_h;/* '<S14>/WGS84 Gravity S-Function' */
  real_T WGS84GravitySFunction_gamma_phi;/* '<S14>/WGS84 Gravity S-Function' */
  real_T IC_FirstOutputTime;           /* '<S13>/IC' */
  real_T IC1_FirstOutputTime;          /* '<S13>/IC1' */
  real_T IC2_FirstOutputTime_d;        /* '<S13>/IC2' */
  real_T IC3_FirstOutputTime_o;        /* '<S13>/IC3' */
  real_T IC4_FirstOutputTime_m;        /* '<S13>/IC4' */
  real_T IC5_FirstOutputTime_d;        /* '<S13>/IC5' */
  real_T otime_PreviousInput;          /* '<S112>/otime' */
  real_T olon_PreviousInput;           /* '<S111>/olon' */
  real_T olat_PreviousInput;           /* '<S110>/olat' */
  real_T oalt_PreviousInput;           /* '<S110>/oalt' */
  real_T NextOutput;                   /* '<S65>/White Noise' */
  real_T NextOutput_a[3];              /* '<S66>/White Noise' */
  real_T TimeStampA_p;                 /* '<S498>/Derivative' */
  real_T LastUAtTimeA_g[3];            /* '<S498>/Derivative' */
  real_T TimeStampB_k;                 /* '<S498>/Derivative' */
  real_T LastUAtTimeB_a[3];            /* '<S498>/Derivative' */
  real_T NextOutput_e[3];              /* '<S534>/White Noise' */
  real_T NextOutput_b[3];              /* '<S548>/White Noise' */
  real_T IC1_FirstOutputTime_k;        /* '<S242>/IC1' */
  real_T IC2_FirstOutputTime_g;        /* '<S242>/IC2' */
  real_T IC9_FirstOutputTime;          /* '<S4>/IC9' */
  real_T Product2_DWORK4[9];           /* '<S292>/Product2' */
  real_T IC13_FirstOutputTime;         /* '<S4>/IC13' */
  real_T IC10_FirstOutputTime;         /* '<S4>/IC10' */
  real_T IC6_FirstOutputTime;          /* '<S4>/IC6' */
  real_T IC_FirstOutputTime_n;         /* '<S4>/IC' */
  real_T IC1_FirstOutputTime_i;        /* '<S4>/IC1' */
  struct {
    void *LoggedData[2];
  } Alphaandbetavisualization_PWORK;   /* '<S7>/Alpha and beta visualization' */

  struct {
    void *LoggedData;
  } Altitude_PWORK;                    /* '<S7>/Altitude' */

  struct {
    void *LoggedData[3];
  } Earthposition_PWORK;               /* '<S7>/Earth position' */

  struct {
    void *LoggedData[3];
  } Euleranglesvisualization_PWORK;    /* '<S7>/Euler angles visualization' */

  struct {
    void *LoggedData[2];
  } Scope_PWORK;                       /* '<S7>/Scope' */

  void *VRSink_PWORK[7];               /* '<S7>/VR Sink' */
  struct {
    void *LoggedData[3];
  } V_bodyvisualization_PWORK;         /* '<S7>/V_body visualization' */

  struct {
    void *LoggedData[3];
  } V_nedvisualization_PWORK;          /* '<S7>/V_ned visualization' */

  struct {
    void *TimePtr;
    void *DataPtr;
    void *RSimInfoPtr;
  } FromWs_PWORK;                      /* '<S463>/FromWs' */

  void *JoystickInput_PWORK[3];        /* '<S5>/Joystick Input' */
  struct {
    void *LoggedData;
  } Scope_PWORK_d;                     /* '<S242>/Scope' */

  int32_T PreLookUpIndexSearch_DWORK1; /* '<S336>/PreLook-Up Index Search' */
  int32_T preAlpha1_DWORK1;            /* '<S242>/preAlpha1' */
  int32_T preBeta_DWORK1;              /* '<S242>/preBeta' */
  int32_T hlookup1_DWORK1;             /* '<S242>/h  look-up1' */
  int32_T lookup1_DWORK1;              /* '<S242>/ look-up1' */
  int32_T p1_DWORK1;                   /* '<S242>/p1' */
  int32_T lookup1_DWORK1_l;            /* '<S242>/look-up1' */
  int32_T frlookup1_DWORK1;            /* '<S242>/fr look-up1' */
  int32_T prelookups1_DWORK1;          /* '<S242>/pre look-ups1' */
  uint32_T PreLookUpIndexSearchaltitude_DWORK1;/* '<S73>/PreLook-Up Index Search  (altitude)' */
  uint32_T PreLookUpIndexSearchprobofexceed_DWORK1;/* '<S73>/PreLook-Up Index Search  (prob of exceed)' */
  uint32_T RandSeed;                   /* '<S65>/White Noise' */
  uint32_T RandSeed_n[3];              /* '<S66>/White Noise' */
  uint32_T RandSeed_g[3];              /* '<S534>/White Noise' */
  uint32_T RandSeed_j[3];              /* '<S548>/White Noise' */
  struct {
    int_T PrevIndex;
  } FromWs_IWORK;                      /* '<S463>/FromWs' */

  int_T JoystickInput_IWORK[5];        /* '<S5>/Joystick Input' */
  int8_T If_ActiveSubsystem;           /* '<S553>/If' */
  int8_T ifHeightMaxlowaltitudeelseifHeightMinisotropicaltitude_ActiveSubsystem;/* '<S60>/if Height < Max low altitude  elseif Height > Min isotropic altitude ' */
  int8_T
    ifHeightMaxlowaltitudeelseifHeightMinisotropicaltitude_ActiveSubsystem_d;/* '<S61>/if Height < Max low altitude  elseif Height > Min isotropic altitude ' */
  uint8_T ForIterator_IterationMarker[7];/* '<S114>/For Iterator' */
  boolean_T IC_FirstOutputTime_m;      /* '<S163>/IC' */
  boolean_T IC1_FirstOutputTime_c;     /* '<S163>/IC1' */
  boolean_T IC2_FirstOutputTime_do;    /* '<S163>/IC2' */
  boolean_T IC3_FirstOutputTime_m;     /* '<S163>/IC3' */
  boolean_T IC4_FirstOutputTime_c;     /* '<S163>/IC4' */
  boolean_T Convertfromgeodetictosphericalcoordinates_MODE;/* '<S104>/Convert from geodetic to  spherical coordinates ' */
  boolean_T Convertfromgeodetictosphericalcoordinates_MODE_n;/* '<S104>/Convert from geodetic to  spherical coordinates' */
  boolean_T Distanceintogustx_MODE;    /* '<S47>/Distance into gust (x)' */
  boolean_T Hpgw_MODE;                 /* '<S55>/Hpgw' */
  boolean_T Hwgwz_MODE;                /* '<S56>/Hwgw(z)' */
  boolean_T Hqgw_MODE;                 /* '<S55>/Hqgw' */
  boolean_T Hvgwz_MODE;                /* '<S56>/Hvgw(z)' */
  boolean_T Hrgw_MODE;                 /* '<S55>/Hrgw' */
  boolean_T Hugwz_MODE;                /* '<S56>/Hugw(z)' */
  boolean_T Quadcopter_MODE;           /* '<S3>/Quadcopter' */
  boolean_T QuadcopterFixedWing_MODE;  /* '<S3>/Quadcopter --> Fixed-Wing' */
  boolean_T FixedWingClimb_MODE;       /* '<S3>/Fixed-Wing Climb' */
  boolean_T FixedWingCruise_MODE;      /* '<S3>/Fixed-Wing - Cruise' */
  boolean_T TakeOff_MODE;              /* '<S3>/TakeOff' */
  boolean_T SpecialcaseNorthSouthGeographicPole_MODE;/* '<S115>/Special case - North//South Geographic Pole' */
  DW_Distanceintogusty_mainV03_56_T Distanceintogustz;/* '<S47>/Distance into gust (z)' */
  DW_Distanceintogusty_mainV03_56_T Distanceintogusty;/* '<S47>/Distance into gust (y)' */
} DW_mainV03_56_T;

/* Continuous states (auto storage) */
typedef struct {
  real_T xeyeze_CSTATE[3];             /* '<S243>/xe,ye,ze' */
  real_T phithetapsi_CSTATE[3];        /* '<S291>/phi theta psi' */
  real_T ubvbwb_CSTATE[3];             /* '<S243>/ub,vb,wb' */
  real_T pqr_CSTATE[3];                /* '<S243>/p,q,r ' */
  real_T TransferFcn_CSTATE;           /* '<S250>/Transfer Fcn' */
  real_T TransferFcn1_CSTATE;          /* '<S250>/Transfer Fcn1' */
  real_T TransferFcn2_CSTATE;          /* '<S250>/Transfer Fcn2' */
  real_T TransferFcn3_CSTATE;          /* '<S250>/Transfer Fcn3' */
  real_T TransferFcn4_CSTATE;          /* '<S250>/Transfer Fcn4' */
  real_T TransferFcn5_CSTATE;          /* '<S250>/Transfer Fcn5' */
  real_T TransferFcn1_CSTATE_k;        /* '<S249>/Transfer Fcn1' */
  real_T TransferFcn4_CSTATE_l;        /* '<S249>/Transfer Fcn4' */
  real_T TransferFcn5_CSTATE_f;        /* '<S249>/Transfer Fcn5' */
  real_T Integrator_CSTATE;            /* '<S247>/Integrator' */
  real_T TransferFcnX_CSTATE[2];       /* '<S538>/Transfer Fcn X' */
  real_T TransferFcnY_CSTATE[2];       /* '<S538>/Transfer Fcn Y' */
  real_T TransferFcnZ_CSTATE[2];       /* '<S538>/Transfer Fcn Z' */
  real_T TransferFcnX_CSTATE_e[2];     /* '<S550>/Transfer Fcn X' */
  real_T TransferFcnY_CSTATE_j[2];     /* '<S550>/Transfer Fcn Y' */
  real_T TransferFcnZ_CSTATE_l[2];     /* '<S550>/Transfer Fcn Z' */
  real_T phithetapsi_CSTATE_e[3];      /* '<S525>/phi theta psi' */
  real_T Integrator_CSTATE_p;          /* '<S233>/Integrator' */
  real_T Filter_CSTATE;                /* '<S233>/Filter' */
  real_T Integrator_CSTATE_i;          /* '<S235>/Integrator' */
  real_T Filter_CSTATE_e;              /* '<S235>/Filter' */
  real_T Integrator_CSTATE_n;          /* '<S234>/Integrator' */
  real_T Filter_CSTATE_b;              /* '<S234>/Filter' */
  real_T Integrator_CSTATE_d;          /* '<S237>/Integrator' */
  real_T Filter_CSTATE_g;              /* '<S237>/Filter' */
  real_T Integrator_CSTATE_k;          /* '<S236>/Integrator' */
  real_T Filter_CSTATE_gp;             /* '<S236>/Filter' */
  real_T Integrator_CSTATE_m;          /* '<S239>/Integrator' */
  real_T Filter_CSTATE_k;              /* '<S239>/Filter' */
  real_T Integrator_CSTATE_py;         /* '<S238>/Integrator' */
  real_T Filter_CSTATE_n;              /* '<S238>/Filter' */
  real_T Integrator_CSTATE_n1;         /* '<S201>/Integrator' */
  real_T Filter_CSTATE_nt;             /* '<S201>/Filter' */
  real_T Integrator_CSTATE_e;          /* '<S207>/Integrator' */
  real_T Filter_CSTATE_m;              /* '<S207>/Filter' */
  real_T Integrator_CSTATE_il;         /* '<S206>/Integrator' */
  real_T Filter_CSTATE_i;              /* '<S206>/Filter' */
  real_T Integrator_CSTATE_mj;         /* '<S209>/Integrator' */
  real_T Filter_CSTATE_kd;             /* '<S209>/Filter' */
  real_T Integrator_CSTATE_dk;         /* '<S208>/Integrator' */
  real_T Filter_CSTATE_ii;             /* '<S208>/Filter' */
  real_T Integrator_CSTATE_l;          /* '<S211>/Integrator' */
  real_T Filter_CSTATE_c;              /* '<S211>/Filter' */
  real_T Integrator_CSTATE_p0;         /* '<S210>/Integrator' */
  real_T Filter_CSTATE_a;              /* '<S210>/Filter' */
  real_T Filter_CSTATE_mr;             /* '<S200>/Filter' */
  real_T Integrator_CSTATE_mx;         /* '<S200>/Integrator' */
  real_T Integrator_CSTATE_g;          /* '<S212>/Integrator' */
  real_T Filter_CSTATE_o;              /* '<S212>/Filter' */
  real_T Integrator_CSTATE_o;          /* '<S202>/Integrator' */
  real_T Filter_CSTATE_mk;             /* '<S202>/Filter' */
  real_T Integrator_CSTATE_io;         /* '<S213>/Integrator' */
  real_T Filter_CSTATE_p;              /* '<S213>/Filter' */
  real_T Integrator_CSTATE_m3;         /* '<S205>/Integrator' */
  real_T Filter_CSTATE_kdq;            /* '<S205>/Filter' */
  real_T Integrator_CSTATE_gq;         /* '<S219>/Integrator' */
  real_T Filter_CSTATE_az;             /* '<S219>/Filter' */
  real_T Integrator_CSTATE_nb;         /* '<S215>/Integrator' */
  real_T Filter_CSTATE_f;              /* '<S215>/Filter' */
  real_T Integrator_CSTATE_j;          /* '<S220>/Integrator' */
  real_T Filter_CSTATE_ok;             /* '<S220>/Filter' */
  real_T Integrator_CSTATE_h;          /* '<S216>/Integrator' */
  real_T Filter_CSTATE_or;             /* '<S216>/Filter' */
  real_T Integrator_CSTATE_b;          /* '<S204>/Integrator' */
  real_T Filter_CSTATE_j;              /* '<S204>/Filter' */
  real_T Integrator_CSTATE_e5;         /* '<S218>/Integrator' */
  real_T Filter_CSTATE_mkr;            /* '<S218>/Filter' */
  real_T Integrator_CSTATE_ni;         /* '<S214>/Integrator' */
  real_T Filter_CSTATE_l;              /* '<S214>/Filter' */
  real_T Integrator_CSTATE_c;          /* '<S217>/Integrator' */
  real_T Filter_CSTATE_j2;             /* '<S217>/Filter' */
  real_T Filter_CSTATE_br;             /* '<S203>/Filter' */
  real_T Integrator_CSTATE_a;          /* '<S203>/Integrator' */
  real_T Filter_CSTATE_ci;             /* '<S179>/Filter' */
  real_T Integrator_CSTATE_nu;         /* '<S179>/Integrator' */
  real_T Integrator_CSTATE_aw;         /* '<S184>/Integrator' */
  real_T Filter_CSTATE_aw;             /* '<S184>/Filter' */
  real_T Integrator_CSTATE_ij;         /* '<S180>/Integrator' */
  real_T Filter_CSTATE_jg;             /* '<S180>/Filter' */
  real_T Integrator_CSTATE_hz;         /* '<S185>/Integrator' */
  real_T Filter_CSTATE_l3;             /* '<S185>/Filter' */
  real_T Integrator_CSTATE_pc;         /* '<S183>/Integrator' */
  real_T Filter_CSTATE_c0;             /* '<S183>/Filter' */
  real_T Integrator_CSTATE_l2;         /* '<S191>/Integrator' */
  real_T Filter_CSTATE_a5;             /* '<S191>/Filter' */
  real_T Integrator_CSTATE_am;         /* '<S187>/Integrator' */
  real_T Filter_CSTATE_mg;             /* '<S187>/Filter' */
  real_T Integrator_CSTATE_dkb;        /* '<S192>/Integrator' */
  real_T Filter_CSTATE_oz;             /* '<S192>/Filter' */
  real_T Integrator_CSTATE_i5;         /* '<S188>/Integrator' */
  real_T Filter_CSTATE_px;             /* '<S188>/Filter' */
  real_T Integrator_CSTATE_ng;         /* '<S182>/Integrator' */
  real_T Filter_CSTATE_d;              /* '<S182>/Filter' */
  real_T Integrator_CSTATE_cp;         /* '<S190>/Integrator' */
  real_T Filter_CSTATE_gi;             /* '<S190>/Filter' */
  real_T Integrator_CSTATE_e1;         /* '<S186>/Integrator' */
  real_T Filter_CSTATE_pc;             /* '<S186>/Filter' */
  real_T Integrator_CSTATE_gw;         /* '<S189>/Integrator' */
  real_T Filter_CSTATE_py;             /* '<S189>/Filter' */
  real_T Filter_CSTATE_ib;             /* '<S181>/Filter' */
  real_T Integrator_CSTATE_l4;         /* '<S181>/Integrator' */
  real_T Integrator_CSTATE_dg;         /* '<S176>/Integrator' */
  real_T Filter_CSTATE_k1;             /* '<S176>/Filter' */
  real_T Integrator_CSTATE_ke;         /* '<S169>/Integrator' */
  real_T Filter_CSTATE_ba;             /* '<S169>/Filter' */
  real_T Integrator_CSTATE_ih;         /* '<S171>/Integrator' */
  real_T Filter_CSTATE_lp;             /* '<S171>/Filter' */
  real_T Integrator_CSTATE_nm;         /* '<S170>/Integrator' */
  real_T Filter_CSTATE_dt;             /* '<S170>/Filter' */
  real_T Integrator_CSTATE_og;         /* '<S173>/Integrator' */
  real_T Filter_CSTATE_iv;             /* '<S173>/Filter' */
  real_T Integrator_CSTATE_ht;         /* '<S172>/Integrator' */
  real_T Filter_CSTATE_k5;             /* '<S172>/Filter' */
  real_T Integrator_CSTATE_cw;         /* '<S175>/Integrator' */
  real_T Filter_CSTATE_lo;             /* '<S175>/Filter' */
  real_T Integrator_CSTATE_ki;         /* '<S174>/Integrator' */
  real_T Filter_CSTATE_kc;             /* '<S174>/Filter' */
  X_Distanceintogusty_mainV03_56_T Distanceintogustz;/* '<S47>/Distance into gust (y)' */
  X_Distanceintogusty_mainV03_56_T Distanceintogusty;/* '<S47>/Distance into gust (y)' */
  real_T DistanceintoGustxLimitedtogustlengthd_CSTATE_a;/* '<S50>/Distance into Gust (x) (Limited to gust length d)' */
} X_mainV03_56_T;

/* State derivatives (auto storage) */
typedef struct {
  real_T xeyeze_CSTATE[3];             /* '<S243>/xe,ye,ze' */
  real_T phithetapsi_CSTATE[3];        /* '<S291>/phi theta psi' */
  real_T ubvbwb_CSTATE[3];             /* '<S243>/ub,vb,wb' */
  real_T pqr_CSTATE[3];                /* '<S243>/p,q,r ' */
  real_T TransferFcn_CSTATE;           /* '<S250>/Transfer Fcn' */
  real_T TransferFcn1_CSTATE;          /* '<S250>/Transfer Fcn1' */
  real_T TransferFcn2_CSTATE;          /* '<S250>/Transfer Fcn2' */
  real_T TransferFcn3_CSTATE;          /* '<S250>/Transfer Fcn3' */
  real_T TransferFcn4_CSTATE;          /* '<S250>/Transfer Fcn4' */
  real_T TransferFcn5_CSTATE;          /* '<S250>/Transfer Fcn5' */
  real_T TransferFcn1_CSTATE_k;        /* '<S249>/Transfer Fcn1' */
  real_T TransferFcn4_CSTATE_l;        /* '<S249>/Transfer Fcn4' */
  real_T TransferFcn5_CSTATE_f;        /* '<S249>/Transfer Fcn5' */
  real_T Integrator_CSTATE;            /* '<S247>/Integrator' */
  real_T TransferFcnX_CSTATE[2];       /* '<S538>/Transfer Fcn X' */
  real_T TransferFcnY_CSTATE[2];       /* '<S538>/Transfer Fcn Y' */
  real_T TransferFcnZ_CSTATE[2];       /* '<S538>/Transfer Fcn Z' */
  real_T TransferFcnX_CSTATE_e[2];     /* '<S550>/Transfer Fcn X' */
  real_T TransferFcnY_CSTATE_j[2];     /* '<S550>/Transfer Fcn Y' */
  real_T TransferFcnZ_CSTATE_l[2];     /* '<S550>/Transfer Fcn Z' */
  real_T phithetapsi_CSTATE_e[3];      /* '<S525>/phi theta psi' */
  real_T Integrator_CSTATE_p;          /* '<S233>/Integrator' */
  real_T Filter_CSTATE;                /* '<S233>/Filter' */
  real_T Integrator_CSTATE_i;          /* '<S235>/Integrator' */
  real_T Filter_CSTATE_e;              /* '<S235>/Filter' */
  real_T Integrator_CSTATE_n;          /* '<S234>/Integrator' */
  real_T Filter_CSTATE_b;              /* '<S234>/Filter' */
  real_T Integrator_CSTATE_d;          /* '<S237>/Integrator' */
  real_T Filter_CSTATE_g;              /* '<S237>/Filter' */
  real_T Integrator_CSTATE_k;          /* '<S236>/Integrator' */
  real_T Filter_CSTATE_gp;             /* '<S236>/Filter' */
  real_T Integrator_CSTATE_m;          /* '<S239>/Integrator' */
  real_T Filter_CSTATE_k;              /* '<S239>/Filter' */
  real_T Integrator_CSTATE_py;         /* '<S238>/Integrator' */
  real_T Filter_CSTATE_n;              /* '<S238>/Filter' */
  real_T Integrator_CSTATE_n1;         /* '<S201>/Integrator' */
  real_T Filter_CSTATE_nt;             /* '<S201>/Filter' */
  real_T Integrator_CSTATE_e;          /* '<S207>/Integrator' */
  real_T Filter_CSTATE_m;              /* '<S207>/Filter' */
  real_T Integrator_CSTATE_il;         /* '<S206>/Integrator' */
  real_T Filter_CSTATE_i;              /* '<S206>/Filter' */
  real_T Integrator_CSTATE_mj;         /* '<S209>/Integrator' */
  real_T Filter_CSTATE_kd;             /* '<S209>/Filter' */
  real_T Integrator_CSTATE_dk;         /* '<S208>/Integrator' */
  real_T Filter_CSTATE_ii;             /* '<S208>/Filter' */
  real_T Integrator_CSTATE_l;          /* '<S211>/Integrator' */
  real_T Filter_CSTATE_c;              /* '<S211>/Filter' */
  real_T Integrator_CSTATE_p0;         /* '<S210>/Integrator' */
  real_T Filter_CSTATE_a;              /* '<S210>/Filter' */
  real_T Filter_CSTATE_mr;             /* '<S200>/Filter' */
  real_T Integrator_CSTATE_mx;         /* '<S200>/Integrator' */
  real_T Integrator_CSTATE_g;          /* '<S212>/Integrator' */
  real_T Filter_CSTATE_o;              /* '<S212>/Filter' */
  real_T Integrator_CSTATE_o;          /* '<S202>/Integrator' */
  real_T Filter_CSTATE_mk;             /* '<S202>/Filter' */
  real_T Integrator_CSTATE_io;         /* '<S213>/Integrator' */
  real_T Filter_CSTATE_p;              /* '<S213>/Filter' */
  real_T Integrator_CSTATE_m3;         /* '<S205>/Integrator' */
  real_T Filter_CSTATE_kdq;            /* '<S205>/Filter' */
  real_T Integrator_CSTATE_gq;         /* '<S219>/Integrator' */
  real_T Filter_CSTATE_az;             /* '<S219>/Filter' */
  real_T Integrator_CSTATE_nb;         /* '<S215>/Integrator' */
  real_T Filter_CSTATE_f;              /* '<S215>/Filter' */
  real_T Integrator_CSTATE_j;          /* '<S220>/Integrator' */
  real_T Filter_CSTATE_ok;             /* '<S220>/Filter' */
  real_T Integrator_CSTATE_h;          /* '<S216>/Integrator' */
  real_T Filter_CSTATE_or;             /* '<S216>/Filter' */
  real_T Integrator_CSTATE_b;          /* '<S204>/Integrator' */
  real_T Filter_CSTATE_j;              /* '<S204>/Filter' */
  real_T Integrator_CSTATE_e5;         /* '<S218>/Integrator' */
  real_T Filter_CSTATE_mkr;            /* '<S218>/Filter' */
  real_T Integrator_CSTATE_ni;         /* '<S214>/Integrator' */
  real_T Filter_CSTATE_l;              /* '<S214>/Filter' */
  real_T Integrator_CSTATE_c;          /* '<S217>/Integrator' */
  real_T Filter_CSTATE_j2;             /* '<S217>/Filter' */
  real_T Filter_CSTATE_br;             /* '<S203>/Filter' */
  real_T Integrator_CSTATE_a;          /* '<S203>/Integrator' */
  real_T Filter_CSTATE_ci;             /* '<S179>/Filter' */
  real_T Integrator_CSTATE_nu;         /* '<S179>/Integrator' */
  real_T Integrator_CSTATE_aw;         /* '<S184>/Integrator' */
  real_T Filter_CSTATE_aw;             /* '<S184>/Filter' */
  real_T Integrator_CSTATE_ij;         /* '<S180>/Integrator' */
  real_T Filter_CSTATE_jg;             /* '<S180>/Filter' */
  real_T Integrator_CSTATE_hz;         /* '<S185>/Integrator' */
  real_T Filter_CSTATE_l3;             /* '<S185>/Filter' */
  real_T Integrator_CSTATE_pc;         /* '<S183>/Integrator' */
  real_T Filter_CSTATE_c0;             /* '<S183>/Filter' */
  real_T Integrator_CSTATE_l2;         /* '<S191>/Integrator' */
  real_T Filter_CSTATE_a5;             /* '<S191>/Filter' */
  real_T Integrator_CSTATE_am;         /* '<S187>/Integrator' */
  real_T Filter_CSTATE_mg;             /* '<S187>/Filter' */
  real_T Integrator_CSTATE_dkb;        /* '<S192>/Integrator' */
  real_T Filter_CSTATE_oz;             /* '<S192>/Filter' */
  real_T Integrator_CSTATE_i5;         /* '<S188>/Integrator' */
  real_T Filter_CSTATE_px;             /* '<S188>/Filter' */
  real_T Integrator_CSTATE_ng;         /* '<S182>/Integrator' */
  real_T Filter_CSTATE_d;              /* '<S182>/Filter' */
  real_T Integrator_CSTATE_cp;         /* '<S190>/Integrator' */
  real_T Filter_CSTATE_gi;             /* '<S190>/Filter' */
  real_T Integrator_CSTATE_e1;         /* '<S186>/Integrator' */
  real_T Filter_CSTATE_pc;             /* '<S186>/Filter' */
  real_T Integrator_CSTATE_gw;         /* '<S189>/Integrator' */
  real_T Filter_CSTATE_py;             /* '<S189>/Filter' */
  real_T Filter_CSTATE_ib;             /* '<S181>/Filter' */
  real_T Integrator_CSTATE_l4;         /* '<S181>/Integrator' */
  real_T Integrator_CSTATE_dg;         /* '<S176>/Integrator' */
  real_T Filter_CSTATE_k1;             /* '<S176>/Filter' */
  real_T Integrator_CSTATE_ke;         /* '<S169>/Integrator' */
  real_T Filter_CSTATE_ba;             /* '<S169>/Filter' */
  real_T Integrator_CSTATE_ih;         /* '<S171>/Integrator' */
  real_T Filter_CSTATE_lp;             /* '<S171>/Filter' */
  real_T Integrator_CSTATE_nm;         /* '<S170>/Integrator' */
  real_T Filter_CSTATE_dt;             /* '<S170>/Filter' */
  real_T Integrator_CSTATE_og;         /* '<S173>/Integrator' */
  real_T Filter_CSTATE_iv;             /* '<S173>/Filter' */
  real_T Integrator_CSTATE_ht;         /* '<S172>/Integrator' */
  real_T Filter_CSTATE_k5;             /* '<S172>/Filter' */
  real_T Integrator_CSTATE_cw;         /* '<S175>/Integrator' */
  real_T Filter_CSTATE_lo;             /* '<S175>/Filter' */
  real_T Integrator_CSTATE_ki;         /* '<S174>/Integrator' */
  real_T Filter_CSTATE_kc;             /* '<S174>/Filter' */
  XDot_Distanceintogusty_mainV03_56_T Distanceintogustz;/* '<S47>/Distance into gust (y)' */
  XDot_Distanceintogusty_mainV03_56_T Distanceintogusty;/* '<S47>/Distance into gust (y)' */
  real_T DistanceintoGustxLimitedtogustlengthd_CSTATE_a;/* '<S50>/Distance into Gust (x) (Limited to gust length d)' */
} XDot_mainV03_56_T;

/* State disabled  */
typedef struct {
  boolean_T xeyeze_CSTATE[3];          /* '<S243>/xe,ye,ze' */
  boolean_T phithetapsi_CSTATE[3];     /* '<S291>/phi theta psi' */
  boolean_T ubvbwb_CSTATE[3];          /* '<S243>/ub,vb,wb' */
  boolean_T pqr_CSTATE[3];             /* '<S243>/p,q,r ' */
  boolean_T TransferFcn_CSTATE;        /* '<S250>/Transfer Fcn' */
  boolean_T TransferFcn1_CSTATE;       /* '<S250>/Transfer Fcn1' */
  boolean_T TransferFcn2_CSTATE;       /* '<S250>/Transfer Fcn2' */
  boolean_T TransferFcn3_CSTATE;       /* '<S250>/Transfer Fcn3' */
  boolean_T TransferFcn4_CSTATE;       /* '<S250>/Transfer Fcn4' */
  boolean_T TransferFcn5_CSTATE;       /* '<S250>/Transfer Fcn5' */
  boolean_T TransferFcn1_CSTATE_k;     /* '<S249>/Transfer Fcn1' */
  boolean_T TransferFcn4_CSTATE_l;     /* '<S249>/Transfer Fcn4' */
  boolean_T TransferFcn5_CSTATE_f;     /* '<S249>/Transfer Fcn5' */
  boolean_T Integrator_CSTATE;         /* '<S247>/Integrator' */
  boolean_T TransferFcnX_CSTATE[2];    /* '<S538>/Transfer Fcn X' */
  boolean_T TransferFcnY_CSTATE[2];    /* '<S538>/Transfer Fcn Y' */
  boolean_T TransferFcnZ_CSTATE[2];    /* '<S538>/Transfer Fcn Z' */
  boolean_T TransferFcnX_CSTATE_e[2];  /* '<S550>/Transfer Fcn X' */
  boolean_T TransferFcnY_CSTATE_j[2];  /* '<S550>/Transfer Fcn Y' */
  boolean_T TransferFcnZ_CSTATE_l[2];  /* '<S550>/Transfer Fcn Z' */
  boolean_T phithetapsi_CSTATE_e[3];   /* '<S525>/phi theta psi' */
  boolean_T Integrator_CSTATE_p;       /* '<S233>/Integrator' */
  boolean_T Filter_CSTATE;             /* '<S233>/Filter' */
  boolean_T Integrator_CSTATE_i;       /* '<S235>/Integrator' */
  boolean_T Filter_CSTATE_e;           /* '<S235>/Filter' */
  boolean_T Integrator_CSTATE_n;       /* '<S234>/Integrator' */
  boolean_T Filter_CSTATE_b;           /* '<S234>/Filter' */
  boolean_T Integrator_CSTATE_d;       /* '<S237>/Integrator' */
  boolean_T Filter_CSTATE_g;           /* '<S237>/Filter' */
  boolean_T Integrator_CSTATE_k;       /* '<S236>/Integrator' */
  boolean_T Filter_CSTATE_gp;          /* '<S236>/Filter' */
  boolean_T Integrator_CSTATE_m;       /* '<S239>/Integrator' */
  boolean_T Filter_CSTATE_k;           /* '<S239>/Filter' */
  boolean_T Integrator_CSTATE_py;      /* '<S238>/Integrator' */
  boolean_T Filter_CSTATE_n;           /* '<S238>/Filter' */
  boolean_T Integrator_CSTATE_n1;      /* '<S201>/Integrator' */
  boolean_T Filter_CSTATE_nt;          /* '<S201>/Filter' */
  boolean_T Integrator_CSTATE_e;       /* '<S207>/Integrator' */
  boolean_T Filter_CSTATE_m;           /* '<S207>/Filter' */
  boolean_T Integrator_CSTATE_il;      /* '<S206>/Integrator' */
  boolean_T Filter_CSTATE_i;           /* '<S206>/Filter' */
  boolean_T Integrator_CSTATE_mj;      /* '<S209>/Integrator' */
  boolean_T Filter_CSTATE_kd;          /* '<S209>/Filter' */
  boolean_T Integrator_CSTATE_dk;      /* '<S208>/Integrator' */
  boolean_T Filter_CSTATE_ii;          /* '<S208>/Filter' */
  boolean_T Integrator_CSTATE_l;       /* '<S211>/Integrator' */
  boolean_T Filter_CSTATE_c;           /* '<S211>/Filter' */
  boolean_T Integrator_CSTATE_p0;      /* '<S210>/Integrator' */
  boolean_T Filter_CSTATE_a;           /* '<S210>/Filter' */
  boolean_T Filter_CSTATE_mr;          /* '<S200>/Filter' */
  boolean_T Integrator_CSTATE_mx;      /* '<S200>/Integrator' */
  boolean_T Integrator_CSTATE_g;       /* '<S212>/Integrator' */
  boolean_T Filter_CSTATE_o;           /* '<S212>/Filter' */
  boolean_T Integrator_CSTATE_o;       /* '<S202>/Integrator' */
  boolean_T Filter_CSTATE_mk;          /* '<S202>/Filter' */
  boolean_T Integrator_CSTATE_io;      /* '<S213>/Integrator' */
  boolean_T Filter_CSTATE_p;           /* '<S213>/Filter' */
  boolean_T Integrator_CSTATE_m3;      /* '<S205>/Integrator' */
  boolean_T Filter_CSTATE_kdq;         /* '<S205>/Filter' */
  boolean_T Integrator_CSTATE_gq;      /* '<S219>/Integrator' */
  boolean_T Filter_CSTATE_az;          /* '<S219>/Filter' */
  boolean_T Integrator_CSTATE_nb;      /* '<S215>/Integrator' */
  boolean_T Filter_CSTATE_f;           /* '<S215>/Filter' */
  boolean_T Integrator_CSTATE_j;       /* '<S220>/Integrator' */
  boolean_T Filter_CSTATE_ok;          /* '<S220>/Filter' */
  boolean_T Integrator_CSTATE_h;       /* '<S216>/Integrator' */
  boolean_T Filter_CSTATE_or;          /* '<S216>/Filter' */
  boolean_T Integrator_CSTATE_b;       /* '<S204>/Integrator' */
  boolean_T Filter_CSTATE_j;           /* '<S204>/Filter' */
  boolean_T Integrator_CSTATE_e5;      /* '<S218>/Integrator' */
  boolean_T Filter_CSTATE_mkr;         /* '<S218>/Filter' */
  boolean_T Integrator_CSTATE_ni;      /* '<S214>/Integrator' */
  boolean_T Filter_CSTATE_l;           /* '<S214>/Filter' */
  boolean_T Integrator_CSTATE_c;       /* '<S217>/Integrator' */
  boolean_T Filter_CSTATE_j2;          /* '<S217>/Filter' */
  boolean_T Filter_CSTATE_br;          /* '<S203>/Filter' */
  boolean_T Integrator_CSTATE_a;       /* '<S203>/Integrator' */
  boolean_T Filter_CSTATE_ci;          /* '<S179>/Filter' */
  boolean_T Integrator_CSTATE_nu;      /* '<S179>/Integrator' */
  boolean_T Integrator_CSTATE_aw;      /* '<S184>/Integrator' */
  boolean_T Filter_CSTATE_aw;          /* '<S184>/Filter' */
  boolean_T Integrator_CSTATE_ij;      /* '<S180>/Integrator' */
  boolean_T Filter_CSTATE_jg;          /* '<S180>/Filter' */
  boolean_T Integrator_CSTATE_hz;      /* '<S185>/Integrator' */
  boolean_T Filter_CSTATE_l3;          /* '<S185>/Filter' */
  boolean_T Integrator_CSTATE_pc;      /* '<S183>/Integrator' */
  boolean_T Filter_CSTATE_c0;          /* '<S183>/Filter' */
  boolean_T Integrator_CSTATE_l2;      /* '<S191>/Integrator' */
  boolean_T Filter_CSTATE_a5;          /* '<S191>/Filter' */
  boolean_T Integrator_CSTATE_am;      /* '<S187>/Integrator' */
  boolean_T Filter_CSTATE_mg;          /* '<S187>/Filter' */
  boolean_T Integrator_CSTATE_dkb;     /* '<S192>/Integrator' */
  boolean_T Filter_CSTATE_oz;          /* '<S192>/Filter' */
  boolean_T Integrator_CSTATE_i5;      /* '<S188>/Integrator' */
  boolean_T Filter_CSTATE_px;          /* '<S188>/Filter' */
  boolean_T Integrator_CSTATE_ng;      /* '<S182>/Integrator' */
  boolean_T Filter_CSTATE_d;           /* '<S182>/Filter' */
  boolean_T Integrator_CSTATE_cp;      /* '<S190>/Integrator' */
  boolean_T Filter_CSTATE_gi;          /* '<S190>/Filter' */
  boolean_T Integrator_CSTATE_e1;      /* '<S186>/Integrator' */
  boolean_T Filter_CSTATE_pc;          /* '<S186>/Filter' */
  boolean_T Integrator_CSTATE_gw;      /* '<S189>/Integrator' */
  boolean_T Filter_CSTATE_py;          /* '<S189>/Filter' */
  boolean_T Filter_CSTATE_ib;          /* '<S181>/Filter' */
  boolean_T Integrator_CSTATE_l4;      /* '<S181>/Integrator' */
  boolean_T Integrator_CSTATE_dg;      /* '<S176>/Integrator' */
  boolean_T Filter_CSTATE_k1;          /* '<S176>/Filter' */
  boolean_T Integrator_CSTATE_ke;      /* '<S169>/Integrator' */
  boolean_T Filter_CSTATE_ba;          /* '<S169>/Filter' */
  boolean_T Integrator_CSTATE_ih;      /* '<S171>/Integrator' */
  boolean_T Filter_CSTATE_lp;          /* '<S171>/Filter' */
  boolean_T Integrator_CSTATE_nm;      /* '<S170>/Integrator' */
  boolean_T Filter_CSTATE_dt;          /* '<S170>/Filter' */
  boolean_T Integrator_CSTATE_og;      /* '<S173>/Integrator' */
  boolean_T Filter_CSTATE_iv;          /* '<S173>/Filter' */
  boolean_T Integrator_CSTATE_ht;      /* '<S172>/Integrator' */
  boolean_T Filter_CSTATE_k5;          /* '<S172>/Filter' */
  boolean_T Integrator_CSTATE_cw;      /* '<S175>/Integrator' */
  boolean_T Filter_CSTATE_lo;          /* '<S175>/Filter' */
  boolean_T Integrator_CSTATE_ki;      /* '<S174>/Integrator' */
  boolean_T Filter_CSTATE_kc;          /* '<S174>/Filter' */
  XDis_Distanceintogusty_mainV03_56_T Distanceintogustz;/* '<S47>/Distance into gust (y)' */
  XDis_Distanceintogusty_mainV03_56_T Distanceintogusty;/* '<S47>/Distance into gust (y)' */
  boolean_T DistanceintoGustxLimitedtogustlengthd_CSTATE_a;/* '<S50>/Distance into Gust (x) (Limited to gust length d)' */
} XDis_mainV03_56_T;

/* Zero-crossing (trigger) state */
typedef struct {
  ZCE_JKFlipFlop_mainV03_56_T JKFlipFlop_jr;/* '<S476>/J-K Flip-Flop' */
  ZCE_JKFlipFlop_mainV03_56_T JKFlipFlop_j;/* '<S475>/J-K Flip-Flop' */
  ZCE_JKFlipFlop_mainV03_56_T JKFlipFlop_f;/* '<S474>/J-K Flip-Flop' */
  ZCE_JKFlipFlop_mainV03_56_T JKFlipFlop_h;/* '<S473>/J-K Flip-Flop' */
  ZCE_JKFlipFlop_mainV03_56_T JKFlipFlop_ej;/* '<S472>/J-K Flip-Flop' */
  ZCE_JKFlipFlop_mainV03_56_T JKFlipFlop_n;/* '<S471>/J-K Flip-Flop' */
  ZCE_JKFlipFlop_mainV03_56_T JKFlipFlop_b;/* '<S470>/J-K Flip-Flop' */
  ZCE_JKFlipFlop_mainV03_56_T JKFlipFlop_p1;/* '<S469>/J-K Flip-Flop' */
  ZCE_JKFlipFlop_mainV03_56_T JKFlipFlop_a;/* '<S468>/J-K Flip-Flop' */
  ZCE_JKFlipFlop_mainV03_56_T JKFlipFlop_p;/* '<S467>/J-K Flip-Flop' */
  ZCE_JKFlipFlop_mainV03_56_T JKFlipFlop_e;/* '<S466>/J-K Flip-Flop' */
  ZCE_JKFlipFlop_mainV03_56_T JKFlipFlop;/* '<S465>/J-K Flip-Flop' */
} PrevZCX_mainV03_56_T;

#ifndef ODE14X_INTG
#define ODE14X_INTG

/* ODE14X Integration Data */
typedef struct {
  real_T *x0;
  real_T *f0;
  real_T *x1start;
  real_T *f1;
  real_T *Delta;
  real_T *E;
  real_T *fac;
  real_T *DFDX;
  real_T *W;
  int_T *pivots;
  real_T *xtmp;
  real_T *ztmp;
  real_T *M;
  real_T *M1;
  real_T *Edot;
  real_T *xdot;
  real_T *fminusMxdot;
  boolean_T isFirstStep;
} ODE14X_IntgData;

#endif

/* Backward compatible GRT Identifiers */
#define rtB                            mainV03_56_B
#define BlockIO                        B_mainV03_56_T
#define rtX                            mainV03_56_X
#define ContinuousStates               X_mainV03_56_T
#define rtXdot                         mainV03_56_XDot
#define StateDerivatives               XDot_mainV03_56_T
#define tXdis                          mainV03_56_XDis
#define StateDisabled                  XDis_mainV03_56_T
#define rtP                            mainV03_56_P
#define Parameters                     P_mainV03_56_T
#define rtDWork                        mainV03_56_DW
#define D_Work                         DW_mainV03_56_T
#define rtPrevZCSigState               mainV03_56_PrevZCX
#define PrevZCSigStates                PrevZCX_mainV03_56_T

/* Parameters for system: '<S47>/Distance into gust (y)' */
struct P_Distanceintogusty_mainV03_56_T_ {
  real_T x_Y0;                         /* Expression: [0]
                                        * Referenced by: '<S51>/x'
                                        */
  real_T DistanceintoGustxLimitedtogustlengthd_IC;/* Expression: 0
                                                   * Referenced by: '<S51>/Distance into Gust (x) (Limited to gust length d) '
                                                   */
  real_T DistanceintoGustxLimitedtogustlengthd_LowerSat;/* Expression: 0
                                                         * Referenced by: '<S51>/Distance into Gust (x) (Limited to gust length d) '
                                                         */
};

/* Parameters for system: '<S465>/J-K Flip-Flop' */
struct P_JKFlipFlop_mainV03_56_T_ {
  boolean_T Q_Y0;                      /* Expression: initial_condition
                                        * Referenced by: '<S477>/Q'
                                        */
  boolean_T Q_Y0_l;                    /* Expression: ~initial_condition
                                        * Referenced by: '<S477>/!Q'
                                        */
};

/* Parameters (auto storage) */
struct P_mainV03_56_T_ {
  struct_PPAnPhCDcFDL0ZPEazAckG groundReaction;/* Variable: groundReaction
                                                * Referenced by:
                                                *   '<S28>/-D1'
                                                *   '<S28>/-K2'
                                                *   '<S29>/-D1'
                                                *   '<S29>/-K2'
                                                *   '<S30>/-D1'
                                                *   '<S30>/-K2'
                                                */
  real_T delta_e;                      /* Variable: delta_e
                                        * Referenced by:
                                        *   '<S164>/Constant'
                                        *   '<S166>/Constant7'
                                        *   '<S168>/Constant'
                                        */
  real_T deltae_degrees;               /* Variable: deltae_degrees
                                        * Referenced by:
                                        *   '<S164>/Constant9'
                                        *   '<S165>/Constant9'
                                        *   '<S168>/Constant9'
                                        */
  real_T deltafl_degrees;              /* Variable: deltafl_degrees
                                        * Referenced by:
                                        *   '<S164>/Constant2'
                                        *   '<S165>/Constant2'
                                        *   '<S168>/Constant2'
                                        */
  real_T deltafr_degrees;              /* Variable: deltafr_degrees
                                        * Referenced by:
                                        *   '<S164>/Constant1'
                                        *   '<S165>/Constant1'
                                        *   '<S168>/Constant1'
                                        */
  real_T deltar_degrees;               /* Variable: deltar_degrees
                                        * Referenced by:
                                        *   '<S164>/Constant10'
                                        *   '<S165>/Constant10'
                                        *   '<S168>/Constant10'
                                        */
  real_T deltat_1;                     /* Variable: deltat_1
                                        * Referenced by:
                                        *   '<S164>/deltat_6'
                                        *   '<S165>/deltat_6'
                                        *   '<S166>/deltat_6'
                                        *   '<S168>/deltat_6'
                                        */
  real_T deltat_2;                     /* Variable: deltat_2
                                        * Referenced by:
                                        *   '<S164>/deltat_7'
                                        *   '<S165>/deltat_7'
                                        *   '<S166>/deltat_7'
                                        *   '<S168>/deltat_7'
                                        */
  real_T deltat_3;                     /* Variable: deltat_3
                                        * Referenced by:
                                        *   '<S164>/deltat_8'
                                        *   '<S165>/deltat_8'
                                        *   '<S166>/deltat_8'
                                        *   '<S168>/deltat_8'
                                        */
  real_T deltat_4;                     /* Variable: deltat_4
                                        * Referenced by:
                                        *   '<S164>/deltat_9'
                                        *   '<S165>/deltat_9'
                                        *   '<S166>/deltat_9'
                                        *   '<S168>/deltat_9'
                                        */
  real_T deltat_5;                     /* Variable: deltat_5
                                        * Referenced by:
                                        *   '<S164>/deltat_10'
                                        *   '<S165>/deltat_10'
                                        *   '<S166>/deltat_10'
                                        *   '<S168>/deltat_10'
                                        */
  real_T PIDController7_D;             /* Mask Parameter: PIDController7_D
                                        * Referenced by: '<S176>/Derivative Gain'
                                        */
  real_T PHeight_D;                    /* Mask Parameter: PHeight_D
                                        * Referenced by: '<S169>/Derivative Gain'
                                        */
  real_T PIDController1_D;             /* Mask Parameter: PIDController1_D
                                        * Referenced by: '<S171>/Derivative Gain'
                                        */
  real_T PIDController_D;              /* Mask Parameter: PIDController_D
                                        * Referenced by: '<S170>/Derivative Gain'
                                        */
  real_T PIDController3_D;             /* Mask Parameter: PIDController3_D
                                        * Referenced by: '<S173>/Derivative Gain'
                                        */
  real_T PIDController2_D;             /* Mask Parameter: PIDController2_D
                                        * Referenced by: '<S172>/Derivative Gain'
                                        */
  real_T PIDController5_D;             /* Mask Parameter: PIDController5_D
                                        * Referenced by: '<S175>/Derivative Gain'
                                        */
  real_T PIDController4_D;             /* Mask Parameter: PIDController4_D
                                        * Referenced by: '<S174>/Derivative Gain'
                                        */
  real_T ESTE2_D;                      /* Mask Parameter: ESTE2_D
                                        * Referenced by: '<S179>/Derivative Gain'
                                        */
  real_T PIDHeightController1_D;       /* Mask Parameter: PIDHeightController1_D
                                        * Referenced by: '<S184>/Derivative Gain'
                                        */
  real_T PIDAltitudeAcceleration_D;    /* Mask Parameter: PIDAltitudeAcceleration_D
                                        * Referenced by: '<S180>/Derivative Gain'
                                        */
  real_T PIDLateralPosition_D;         /* Mask Parameter: PIDLateralPosition_D
                                        * Referenced by: '<S185>/Derivative Gain'
                                        */
  real_T PIDController1_D_o;           /* Mask Parameter: PIDController1_D_o
                                        * Referenced by: '<S183>/Derivative Gain'
                                        */
  real_T ProportionalRoll_D;           /* Mask Parameter: ProportionalRoll_D
                                        * Referenced by: '<S191>/Derivative Gain'
                                        */
  real_T PIDRollRate_D;                /* Mask Parameter: PIDRollRate_D
                                        * Referenced by: '<S187>/Derivative Gain'
                                        */
  real_T ProportionalYaw1_D;           /* Mask Parameter: ProportionalYaw1_D
                                        * Referenced by: '<S192>/Derivative Gain'
                                        */
  real_T PIDYawRate1_D;                /* Mask Parameter: PIDYawRate1_D
                                        * Referenced by: '<S188>/Derivative Gain'
                                        */
  real_T PIDController_D_n;            /* Mask Parameter: PIDController_D_n
                                        * Referenced by: '<S182>/Derivative Gain'
                                        */
  real_T ProportionalPitch2_D;         /* Mask Parameter: ProportionalPitch2_D
                                        * Referenced by: '<S190>/Derivative Gain'
                                        */
  real_T PIDPitchRate2_D;              /* Mask Parameter: PIDPitchRate2_D
                                        * Referenced by: '<S186>/Derivative Gain'
                                        */
  real_T ProportionalAltitude_D;       /* Mask Parameter: ProportionalAltitude_D
                                        * Referenced by: '<S189>/Derivative Gain'
                                        */
  real_T PIDAltitudeRate_D;            /* Mask Parameter: PIDAltitudeRate_D
                                        * Referenced by: '<S181>/Derivative Gain'
                                        */
  real_T PHeight1_D;                   /* Mask Parameter: PHeight1_D
                                        * Referenced by: '<S201>/Derivative Gain'
                                        */
  real_T PIDController5_D_f;           /* Mask Parameter: PIDController5_D_f
                                        * Referenced by: '<S207>/Derivative Gain'
                                        */
  real_T PIDController4_D_n;           /* Mask Parameter: PIDController4_D_n
                                        * Referenced by: '<S206>/Derivative Gain'
                                        */
  real_T PIDController7_D_i;           /* Mask Parameter: PIDController7_D_i
                                        * Referenced by: '<S209>/Derivative Gain'
                                        */
  real_T PIDController6_D;             /* Mask Parameter: PIDController6_D
                                        * Referenced by: '<S208>/Derivative Gain'
                                        */
  real_T PIDController9_D;             /* Mask Parameter: PIDController9_D
                                        * Referenced by: '<S211>/Derivative Gain'
                                        */
  real_T PIDController8_D;             /* Mask Parameter: PIDController8_D
                                        * Referenced by: '<S210>/Derivative Gain'
                                        */
  real_T ESTE2_D_i;                    /* Mask Parameter: ESTE2_D_i
                                        * Referenced by: '<S200>/Derivative Gain'
                                        */
  real_T PIDHeightController1_D_l;     /* Mask Parameter: PIDHeightController1_D_l
                                        * Referenced by: '<S212>/Derivative Gain'
                                        */
  real_T PIDAltitudeAcceleration_D_l;  /* Mask Parameter: PIDAltitudeAcceleration_D_l
                                        * Referenced by: '<S202>/Derivative Gain'
                                        */
  real_T PIDLateralPosition_D_i;       /* Mask Parameter: PIDLateralPosition_D_i
                                        * Referenced by: '<S213>/Derivative Gain'
                                        */
  real_T PIDController1_D_m;           /* Mask Parameter: PIDController1_D_m
                                        * Referenced by: '<S205>/Derivative Gain'
                                        */
  real_T ProportionalRoll_D_e;         /* Mask Parameter: ProportionalRoll_D_e
                                        * Referenced by: '<S219>/Derivative Gain'
                                        */
  real_T PIDRollRate_D_l;              /* Mask Parameter: PIDRollRate_D_l
                                        * Referenced by: '<S215>/Derivative Gain'
                                        */
  real_T ProportionalYaw1_D_h;         /* Mask Parameter: ProportionalYaw1_D_h
                                        * Referenced by: '<S220>/Derivative Gain'
                                        */
  real_T PIDYawRate1_D_e;              /* Mask Parameter: PIDYawRate1_D_e
                                        * Referenced by: '<S216>/Derivative Gain'
                                        */
  real_T PIDController_D_b;            /* Mask Parameter: PIDController_D_b
                                        * Referenced by: '<S204>/Derivative Gain'
                                        */
  real_T ProportionalPitch2_D_f;       /* Mask Parameter: ProportionalPitch2_D_f
                                        * Referenced by: '<S218>/Derivative Gain'
                                        */
  real_T PIDPitchRate2_D_d;            /* Mask Parameter: PIDPitchRate2_D_d
                                        * Referenced by: '<S214>/Derivative Gain'
                                        */
  real_T ProportionalAltitude_D_h;     /* Mask Parameter: ProportionalAltitude_D_h
                                        * Referenced by: '<S217>/Derivative Gain'
                                        */
  real_T PIDAltitudeRate_D_h;          /* Mask Parameter: PIDAltitudeRate_D_h
                                        * Referenced by: '<S203>/Derivative Gain'
                                        */
  real_T PHeight_D_n;                  /* Mask Parameter: PHeight_D_n
                                        * Referenced by: '<S233>/Derivative Gain'
                                        */
  real_T PIDController1_D_a;           /* Mask Parameter: PIDController1_D_a
                                        * Referenced by: '<S235>/Derivative Gain'
                                        */
  real_T PIDController_D_a;            /* Mask Parameter: PIDController_D_a
                                        * Referenced by: '<S234>/Derivative Gain'
                                        */
  real_T PIDController3_D_g;           /* Mask Parameter: PIDController3_D_g
                                        * Referenced by: '<S237>/Derivative Gain'
                                        */
  real_T PIDController2_D_b;           /* Mask Parameter: PIDController2_D_b
                                        * Referenced by: '<S236>/Derivative Gain'
                                        */
  real_T PIDController5_D_e;           /* Mask Parameter: PIDController5_D_e
                                        * Referenced by: '<S239>/Derivative Gain'
                                        */
  real_T PIDController4_D_nv;          /* Mask Parameter: PIDController4_D_nv
                                        * Referenced by: '<S238>/Derivative Gain'
                                        */
  real_T DiscreteWindGustModel_Gx;     /* Mask Parameter: DiscreteWindGustModel_Gx
                                        * Referenced by: '<S47>/Constant'
                                        */
  real_T DiscreteWindGustModel_Gy;     /* Mask Parameter: DiscreteWindGustModel_Gy
                                        * Referenced by: '<S47>/Constant1'
                                        */
  real_T DiscreteWindGustModel_Gz;     /* Mask Parameter: DiscreteWindGustModel_Gz
                                        * Referenced by: '<S47>/Constant2'
                                        */
  real_T PHeight_I;                    /* Mask Parameter: PHeight_I
                                        * Referenced by: '<S169>/Integral Gain'
                                        */
  real_T PIDController_I;              /* Mask Parameter: PIDController_I
                                        * Referenced by: '<S170>/Integral Gain'
                                        */
  real_T PIDController1_I;             /* Mask Parameter: PIDController1_I
                                        * Referenced by: '<S171>/Integral Gain'
                                        */
  real_T PIDController2_I;             /* Mask Parameter: PIDController2_I
                                        * Referenced by: '<S172>/Integral Gain'
                                        */
  real_T PIDController3_I;             /* Mask Parameter: PIDController3_I
                                        * Referenced by: '<S173>/Integral Gain'
                                        */
  real_T PIDController4_I;             /* Mask Parameter: PIDController4_I
                                        * Referenced by: '<S174>/Integral Gain'
                                        */
  real_T PIDController5_I;             /* Mask Parameter: PIDController5_I
                                        * Referenced by: '<S175>/Integral Gain'
                                        */
  real_T PIDController7_I;             /* Mask Parameter: PIDController7_I
                                        * Referenced by: '<S176>/Integral Gain'
                                        */
  real_T ESTE2_I;                      /* Mask Parameter: ESTE2_I
                                        * Referenced by: '<S179>/Integral Gain'
                                        */
  real_T PIDAltitudeAcceleration_I;    /* Mask Parameter: PIDAltitudeAcceleration_I
                                        * Referenced by: '<S180>/Integral Gain'
                                        */
  real_T PIDAltitudeRate_I;            /* Mask Parameter: PIDAltitudeRate_I
                                        * Referenced by: '<S181>/Integral Gain'
                                        */
  real_T PIDController_I_f;            /* Mask Parameter: PIDController_I_f
                                        * Referenced by: '<S182>/Integral Gain'
                                        */
  real_T PIDController1_I_l;           /* Mask Parameter: PIDController1_I_l
                                        * Referenced by: '<S183>/Integral Gain'
                                        */
  real_T PIDHeightController1_I;       /* Mask Parameter: PIDHeightController1_I
                                        * Referenced by: '<S184>/Integral Gain'
                                        */
  real_T PIDLateralPosition_I;         /* Mask Parameter: PIDLateralPosition_I
                                        * Referenced by: '<S185>/Integral Gain'
                                        */
  real_T PIDPitchRate2_I;              /* Mask Parameter: PIDPitchRate2_I
                                        * Referenced by: '<S186>/Integral Gain'
                                        */
  real_T PIDRollRate_I;                /* Mask Parameter: PIDRollRate_I
                                        * Referenced by: '<S187>/Integral Gain'
                                        */
  real_T PIDYawRate1_I;                /* Mask Parameter: PIDYawRate1_I
                                        * Referenced by: '<S188>/Integral Gain'
                                        */
  real_T ProportionalAltitude_I;       /* Mask Parameter: ProportionalAltitude_I
                                        * Referenced by: '<S189>/Integral Gain'
                                        */
  real_T ProportionalPitch2_I;         /* Mask Parameter: ProportionalPitch2_I
                                        * Referenced by: '<S190>/Integral Gain'
                                        */
  real_T ProportionalRoll_I;           /* Mask Parameter: ProportionalRoll_I
                                        * Referenced by: '<S191>/Integral Gain'
                                        */
  real_T ProportionalYaw1_I;           /* Mask Parameter: ProportionalYaw1_I
                                        * Referenced by: '<S192>/Integral Gain'
                                        */
  real_T ESTE2_I_p;                    /* Mask Parameter: ESTE2_I_p
                                        * Referenced by: '<S200>/Integral Gain'
                                        */
  real_T PHeight1_I;                   /* Mask Parameter: PHeight1_I
                                        * Referenced by: '<S201>/Integral Gain'
                                        */
  real_T PIDAltitudeAcceleration_I_d;  /* Mask Parameter: PIDAltitudeAcceleration_I_d
                                        * Referenced by: '<S202>/Integral Gain'
                                        */
  real_T PIDAltitudeRate_I_p;          /* Mask Parameter: PIDAltitudeRate_I_p
                                        * Referenced by: '<S203>/Integral Gain'
                                        */
  real_T PIDController_I_b;            /* Mask Parameter: PIDController_I_b
                                        * Referenced by: '<S204>/Integral Gain'
                                        */
  real_T PIDController1_I_i;           /* Mask Parameter: PIDController1_I_i
                                        * Referenced by: '<S205>/Integral Gain'
                                        */
  real_T PIDController4_I_f;           /* Mask Parameter: PIDController4_I_f
                                        * Referenced by: '<S206>/Integral Gain'
                                        */
  real_T PIDController5_I_p;           /* Mask Parameter: PIDController5_I_p
                                        * Referenced by: '<S207>/Integral Gain'
                                        */
  real_T PIDController6_I;             /* Mask Parameter: PIDController6_I
                                        * Referenced by: '<S208>/Integral Gain'
                                        */
  real_T PIDController7_I_a;           /* Mask Parameter: PIDController7_I_a
                                        * Referenced by: '<S209>/Integral Gain'
                                        */
  real_T PIDController8_I;             /* Mask Parameter: PIDController8_I
                                        * Referenced by: '<S210>/Integral Gain'
                                        */
  real_T PIDController9_I;             /* Mask Parameter: PIDController9_I
                                        * Referenced by: '<S211>/Integral Gain'
                                        */
  real_T PIDHeightController1_I_n;     /* Mask Parameter: PIDHeightController1_I_n
                                        * Referenced by: '<S212>/Integral Gain'
                                        */
  real_T PIDLateralPosition_I_n;       /* Mask Parameter: PIDLateralPosition_I_n
                                        * Referenced by: '<S213>/Integral Gain'
                                        */
  real_T PIDPitchRate2_I_l;            /* Mask Parameter: PIDPitchRate2_I_l
                                        * Referenced by: '<S214>/Integral Gain'
                                        */
  real_T PIDRollRate_I_n;              /* Mask Parameter: PIDRollRate_I_n
                                        * Referenced by: '<S215>/Integral Gain'
                                        */
  real_T PIDYawRate1_I_j;              /* Mask Parameter: PIDYawRate1_I_j
                                        * Referenced by: '<S216>/Integral Gain'
                                        */
  real_T ProportionalAltitude_I_h;     /* Mask Parameter: ProportionalAltitude_I_h
                                        * Referenced by: '<S217>/Integral Gain'
                                        */
  real_T ProportionalPitch2_I_d;       /* Mask Parameter: ProportionalPitch2_I_d
                                        * Referenced by: '<S218>/Integral Gain'
                                        */
  real_T ProportionalRoll_I_b;         /* Mask Parameter: ProportionalRoll_I_b
                                        * Referenced by: '<S219>/Integral Gain'
                                        */
  real_T ProportionalYaw1_I_a;         /* Mask Parameter: ProportionalYaw1_I_a
                                        * Referenced by: '<S220>/Integral Gain'
                                        */
  real_T PHeight_I_f;                  /* Mask Parameter: PHeight_I_f
                                        * Referenced by: '<S233>/Integral Gain'
                                        */
  real_T PIDController_I_c;            /* Mask Parameter: PIDController_I_c
                                        * Referenced by: '<S234>/Integral Gain'
                                        */
  real_T PIDController1_I_n;           /* Mask Parameter: PIDController1_I_n
                                        * Referenced by: '<S235>/Integral Gain'
                                        */
  real_T PIDController2_I_g;           /* Mask Parameter: PIDController2_I_g
                                        * Referenced by: '<S236>/Integral Gain'
                                        */
  real_T PIDController3_I_f;           /* Mask Parameter: PIDController3_I_f
                                        * Referenced by: '<S237>/Integral Gain'
                                        */
  real_T PIDController4_I_o;           /* Mask Parameter: PIDController4_I_o
                                        * Referenced by: '<S238>/Integral Gain'
                                        */
  real_T PIDController5_I_b;           /* Mask Parameter: PIDController5_I_b
                                        * Referenced by: '<S239>/Integral Gain'
                                        */
  real_T FlatEarthtoLLA_LL0[2];        /* Mask Parameter: FlatEarthtoLLA_LL0
                                        * Referenced by: '<S244>/initial_pos'
                                        */
  real_T LLAtoFlatEarth_LL0[2];        /* Mask Parameter: LLAtoFlatEarth_LL0
                                        * Referenced by: '<S502>/initial_pos'
                                        */
  real_T DrydenWindTurbulenceModelDiscreteqr_L_high;/* Mask Parameter: DrydenWindTurbulenceModelDiscreteqr_L_high
                                                     * Referenced by: '<S92>/Medium//High Altitude'
                                                     */
  real_T PIDController7_LowerSaturationLimit;/* Mask Parameter: PIDController7_LowerSaturationLimit
                                              * Referenced by: '<S176>/Saturate'
                                              */
  real_T PHeight_LowerSaturationLimit; /* Mask Parameter: PHeight_LowerSaturationLimit
                                        * Referenced by: '<S169>/Saturate'
                                        */
  real_T PIDController_LowerSaturationLimit;/* Mask Parameter: PIDController_LowerSaturationLimit
                                             * Referenced by: '<S170>/Saturate'
                                             */
  real_T ProportionalAltitude_LowerSaturationLimit;/* Mask Parameter: ProportionalAltitude_LowerSaturationLimit
                                                    * Referenced by: '<S189>/Saturate'
                                                    */
  real_T PHeight1_LowerSaturationLimit;/* Mask Parameter: PHeight1_LowerSaturationLimit
                                        * Referenced by: '<S201>/Saturate'
                                        */
  real_T PIDController4_LowerSaturationLimit;/* Mask Parameter: PIDController4_LowerSaturationLimit
                                              * Referenced by: '<S206>/Saturate'
                                              */
  real_T ProportionalAltitude_LowerSaturationLimit_c;/* Mask Parameter: ProportionalAltitude_LowerSaturationLimit_c
                                                      * Referenced by: '<S217>/Saturate'
                                                      */
  real_T PHeight_LowerSaturationLimit_g;/* Mask Parameter: PHeight_LowerSaturationLimit_g
                                         * Referenced by: '<S233>/Saturate'
                                         */
  real_T PIDController_LowerSaturationLimit_l;/* Mask Parameter: PIDController_LowerSaturationLimit_l
                                               * Referenced by: '<S234>/Saturate'
                                               */
  real_T PIDController7_N;             /* Mask Parameter: PIDController7_N
                                        * Referenced by: '<S176>/Filter Coefficient'
                                        */
  real_T PHeight_N;                    /* Mask Parameter: PHeight_N
                                        * Referenced by: '<S169>/Filter Coefficient'
                                        */
  real_T PIDController1_N;             /* Mask Parameter: PIDController1_N
                                        * Referenced by: '<S171>/Filter Coefficient'
                                        */
  real_T PIDController_N;              /* Mask Parameter: PIDController_N
                                        * Referenced by: '<S170>/Filter Coefficient'
                                        */
  real_T PIDController3_N;             /* Mask Parameter: PIDController3_N
                                        * Referenced by: '<S173>/Filter Coefficient'
                                        */
  real_T PIDController2_N;             /* Mask Parameter: PIDController2_N
                                        * Referenced by: '<S172>/Filter Coefficient'
                                        */
  real_T PIDController5_N;             /* Mask Parameter: PIDController5_N
                                        * Referenced by: '<S175>/Filter Coefficient'
                                        */
  real_T PIDController4_N;             /* Mask Parameter: PIDController4_N
                                        * Referenced by: '<S174>/Filter Coefficient'
                                        */
  real_T ESTE2_N;                      /* Mask Parameter: ESTE2_N
                                        * Referenced by: '<S179>/Filter Coefficient'
                                        */
  real_T PIDHeightController1_N;       /* Mask Parameter: PIDHeightController1_N
                                        * Referenced by: '<S184>/Filter Coefficient'
                                        */
  real_T PIDAltitudeAcceleration_N;    /* Mask Parameter: PIDAltitudeAcceleration_N
                                        * Referenced by: '<S180>/Filter Coefficient'
                                        */
  real_T PIDLateralPosition_N;         /* Mask Parameter: PIDLateralPosition_N
                                        * Referenced by: '<S185>/Filter Coefficient'
                                        */
  real_T PIDController1_N_p;           /* Mask Parameter: PIDController1_N_p
                                        * Referenced by: '<S183>/Filter Coefficient'
                                        */
  real_T ProportionalRoll_N;           /* Mask Parameter: ProportionalRoll_N
                                        * Referenced by: '<S191>/Filter Coefficient'
                                        */
  real_T PIDRollRate_N;                /* Mask Parameter: PIDRollRate_N
                                        * Referenced by: '<S187>/Filter Coefficient'
                                        */
  real_T ProportionalYaw1_N;           /* Mask Parameter: ProportionalYaw1_N
                                        * Referenced by: '<S192>/Filter Coefficient'
                                        */
  real_T PIDYawRate1_N;                /* Mask Parameter: PIDYawRate1_N
                                        * Referenced by: '<S188>/Filter Coefficient'
                                        */
  real_T PIDController_N_d;            /* Mask Parameter: PIDController_N_d
                                        * Referenced by: '<S182>/Filter Coefficient'
                                        */
  real_T ProportionalPitch2_N;         /* Mask Parameter: ProportionalPitch2_N
                                        * Referenced by: '<S190>/Filter Coefficient'
                                        */
  real_T PIDPitchRate2_N;              /* Mask Parameter: PIDPitchRate2_N
                                        * Referenced by: '<S186>/Filter Coefficient'
                                        */
  real_T ProportionalAltitude_N;       /* Mask Parameter: ProportionalAltitude_N
                                        * Referenced by: '<S189>/Filter Coefficient'
                                        */
  real_T PIDAltitudeRate_N;            /* Mask Parameter: PIDAltitudeRate_N
                                        * Referenced by: '<S181>/Filter Coefficient'
                                        */
  real_T PHeight1_N;                   /* Mask Parameter: PHeight1_N
                                        * Referenced by: '<S201>/Filter Coefficient'
                                        */
  real_T PIDController5_N_p;           /* Mask Parameter: PIDController5_N_p
                                        * Referenced by: '<S207>/Filter Coefficient'
                                        */
  real_T PIDController4_N_p;           /* Mask Parameter: PIDController4_N_p
                                        * Referenced by: '<S206>/Filter Coefficient'
                                        */
  real_T PIDController7_N_j;           /* Mask Parameter: PIDController7_N_j
                                        * Referenced by: '<S209>/Filter Coefficient'
                                        */
  real_T PIDController6_N;             /* Mask Parameter: PIDController6_N
                                        * Referenced by: '<S208>/Filter Coefficient'
                                        */
  real_T PIDController9_N;             /* Mask Parameter: PIDController9_N
                                        * Referenced by: '<S211>/Filter Coefficient'
                                        */
  real_T PIDController8_N;             /* Mask Parameter: PIDController8_N
                                        * Referenced by: '<S210>/Filter Coefficient'
                                        */
  real_T ESTE2_N_d;                    /* Mask Parameter: ESTE2_N_d
                                        * Referenced by: '<S200>/Filter Coefficient'
                                        */
  real_T PIDHeightController1_N_l;     /* Mask Parameter: PIDHeightController1_N_l
                                        * Referenced by: '<S212>/Filter Coefficient'
                                        */
  real_T PIDAltitudeAcceleration_N_l;  /* Mask Parameter: PIDAltitudeAcceleration_N_l
                                        * Referenced by: '<S202>/Filter Coefficient'
                                        */
  real_T PIDLateralPosition_N_p;       /* Mask Parameter: PIDLateralPosition_N_p
                                        * Referenced by: '<S213>/Filter Coefficient'
                                        */
  real_T PIDController1_N_c;           /* Mask Parameter: PIDController1_N_c
                                        * Referenced by: '<S205>/Filter Coefficient'
                                        */
  real_T ProportionalRoll_N_k;         /* Mask Parameter: ProportionalRoll_N_k
                                        * Referenced by: '<S219>/Filter Coefficient'
                                        */
  real_T PIDRollRate_N_d;              /* Mask Parameter: PIDRollRate_N_d
                                        * Referenced by: '<S215>/Filter Coefficient'
                                        */
  real_T ProportionalYaw1_N_j;         /* Mask Parameter: ProportionalYaw1_N_j
                                        * Referenced by: '<S220>/Filter Coefficient'
                                        */
  real_T PIDYawRate1_N_h;              /* Mask Parameter: PIDYawRate1_N_h
                                        * Referenced by: '<S216>/Filter Coefficient'
                                        */
  real_T PIDController_N_g;            /* Mask Parameter: PIDController_N_g
                                        * Referenced by: '<S204>/Filter Coefficient'
                                        */
  real_T ProportionalPitch2_N_b;       /* Mask Parameter: ProportionalPitch2_N_b
                                        * Referenced by: '<S218>/Filter Coefficient'
                                        */
  real_T PIDPitchRate2_N_j;            /* Mask Parameter: PIDPitchRate2_N_j
                                        * Referenced by: '<S214>/Filter Coefficient'
                                        */
  real_T ProportionalAltitude_N_e;     /* Mask Parameter: ProportionalAltitude_N_e
                                        * Referenced by: '<S217>/Filter Coefficient'
                                        */
  real_T PIDAltitudeRate_N_k;          /* Mask Parameter: PIDAltitudeRate_N_k
                                        * Referenced by: '<S203>/Filter Coefficient'
                                        */
  real_T PHeight_N_h;                  /* Mask Parameter: PHeight_N_h
                                        * Referenced by: '<S233>/Filter Coefficient'
                                        */
  real_T PIDController1_N_m;           /* Mask Parameter: PIDController1_N_m
                                        * Referenced by: '<S235>/Filter Coefficient'
                                        */
  real_T PIDController_N_a;            /* Mask Parameter: PIDController_N_a
                                        * Referenced by: '<S234>/Filter Coefficient'
                                        */
  real_T PIDController3_N_h;           /* Mask Parameter: PIDController3_N_h
                                        * Referenced by: '<S237>/Filter Coefficient'
                                        */
  real_T PIDController2_N_m;           /* Mask Parameter: PIDController2_N_m
                                        * Referenced by: '<S236>/Filter Coefficient'
                                        */
  real_T PIDController5_N_m;           /* Mask Parameter: PIDController5_N_m
                                        * Referenced by: '<S239>/Filter Coefficient'
                                        */
  real_T PIDController4_N_f;           /* Mask Parameter: PIDController4_N_f
                                        * Referenced by: '<S238>/Filter Coefficient'
                                        */
  real_T PIDController7_P;             /* Mask Parameter: PIDController7_P
                                        * Referenced by: '<S176>/Proportional Gain'
                                        */
  real_T PHeight_P;                    /* Mask Parameter: PHeight_P
                                        * Referenced by: '<S169>/Proportional Gain'
                                        */
  real_T PIDController1_P;             /* Mask Parameter: PIDController1_P
                                        * Referenced by: '<S171>/Proportional Gain'
                                        */
  real_T PIDController_P;              /* Mask Parameter: PIDController_P
                                        * Referenced by: '<S170>/Proportional Gain'
                                        */
  real_T PIDController3_P;             /* Mask Parameter: PIDController3_P
                                        * Referenced by: '<S173>/Proportional Gain'
                                        */
  real_T PIDController2_P;             /* Mask Parameter: PIDController2_P
                                        * Referenced by: '<S172>/Proportional Gain'
                                        */
  real_T PIDController5_P;             /* Mask Parameter: PIDController5_P
                                        * Referenced by: '<S175>/Proportional Gain'
                                        */
  real_T PIDController4_P;             /* Mask Parameter: PIDController4_P
                                        * Referenced by: '<S174>/Proportional Gain'
                                        */
  real_T ESTE2_P;                      /* Mask Parameter: ESTE2_P
                                        * Referenced by: '<S179>/Proportional Gain'
                                        */
  real_T PIDHeightController1_P;       /* Mask Parameter: PIDHeightController1_P
                                        * Referenced by: '<S184>/Proportional Gain'
                                        */
  real_T PIDAltitudeAcceleration_P;    /* Mask Parameter: PIDAltitudeAcceleration_P
                                        * Referenced by: '<S180>/Proportional Gain'
                                        */
  real_T PIDLateralPosition_P;         /* Mask Parameter: PIDLateralPosition_P
                                        * Referenced by: '<S185>/Proportional Gain'
                                        */
  real_T PIDController1_P_h;           /* Mask Parameter: PIDController1_P_h
                                        * Referenced by: '<S183>/Proportional Gain'
                                        */
  real_T ProportionalRoll_P;           /* Mask Parameter: ProportionalRoll_P
                                        * Referenced by: '<S191>/Proportional Gain'
                                        */
  real_T PIDRollRate_P;                /* Mask Parameter: PIDRollRate_P
                                        * Referenced by: '<S187>/Proportional Gain'
                                        */
  real_T ProportionalYaw1_P;           /* Mask Parameter: ProportionalYaw1_P
                                        * Referenced by: '<S192>/Proportional Gain'
                                        */
  real_T PIDYawRate1_P;                /* Mask Parameter: PIDYawRate1_P
                                        * Referenced by: '<S188>/Proportional Gain'
                                        */
  real_T PIDController_P_n;            /* Mask Parameter: PIDController_P_n
                                        * Referenced by: '<S182>/Proportional Gain'
                                        */
  real_T ProportionalPitch2_P;         /* Mask Parameter: ProportionalPitch2_P
                                        * Referenced by: '<S190>/Proportional Gain'
                                        */
  real_T PIDPitchRate2_P;              /* Mask Parameter: PIDPitchRate2_P
                                        * Referenced by: '<S186>/Proportional Gain'
                                        */
  real_T ProportionalAltitude_P;       /* Mask Parameter: ProportionalAltitude_P
                                        * Referenced by: '<S189>/Proportional Gain'
                                        */
  real_T PHeight1_P;                   /* Mask Parameter: PHeight1_P
                                        * Referenced by: '<S201>/Proportional Gain'
                                        */
  real_T PIDController5_P_n;           /* Mask Parameter: PIDController5_P_n
                                        * Referenced by: '<S207>/Proportional Gain'
                                        */
  real_T PIDController4_P_k;           /* Mask Parameter: PIDController4_P_k
                                        * Referenced by: '<S206>/Proportional Gain'
                                        */
  real_T PIDController7_P_a;           /* Mask Parameter: PIDController7_P_a
                                        * Referenced by: '<S209>/Proportional Gain'
                                        */
  real_T PIDController6_P;             /* Mask Parameter: PIDController6_P
                                        * Referenced by: '<S208>/Proportional Gain'
                                        */
  real_T PIDController9_P;             /* Mask Parameter: PIDController9_P
                                        * Referenced by: '<S211>/Proportional Gain'
                                        */
  real_T PIDController8_P;             /* Mask Parameter: PIDController8_P
                                        * Referenced by: '<S210>/Proportional Gain'
                                        */
  real_T ESTE2_P_p;                    /* Mask Parameter: ESTE2_P_p
                                        * Referenced by: '<S200>/Proportional Gain'
                                        */
  real_T PIDHeightController1_P_h;     /* Mask Parameter: PIDHeightController1_P_h
                                        * Referenced by: '<S212>/Proportional Gain'
                                        */
  real_T PIDAltitudeAcceleration_P_f;  /* Mask Parameter: PIDAltitudeAcceleration_P_f
                                        * Referenced by: '<S202>/Proportional Gain'
                                        */
  real_T PIDLateralPosition_P_g;       /* Mask Parameter: PIDLateralPosition_P_g
                                        * Referenced by: '<S213>/Proportional Gain'
                                        */
  real_T PIDController1_P_c;           /* Mask Parameter: PIDController1_P_c
                                        * Referenced by: '<S205>/Proportional Gain'
                                        */
  real_T ProportionalRoll_P_g;         /* Mask Parameter: ProportionalRoll_P_g
                                        * Referenced by: '<S219>/Proportional Gain'
                                        */
  real_T PIDRollRate_P_j;              /* Mask Parameter: PIDRollRate_P_j
                                        * Referenced by: '<S215>/Proportional Gain'
                                        */
  real_T ProportionalYaw1_P_p;         /* Mask Parameter: ProportionalYaw1_P_p
                                        * Referenced by: '<S220>/Proportional Gain'
                                        */
  real_T PIDYawRate1_P_g;              /* Mask Parameter: PIDYawRate1_P_g
                                        * Referenced by: '<S216>/Proportional Gain'
                                        */
  real_T PIDController_P_g;            /* Mask Parameter: PIDController_P_g
                                        * Referenced by: '<S204>/Proportional Gain'
                                        */
  real_T ProportionalPitch2_P_i;       /* Mask Parameter: ProportionalPitch2_P_i
                                        * Referenced by: '<S218>/Proportional Gain'
                                        */
  real_T PIDPitchRate2_P_h;            /* Mask Parameter: PIDPitchRate2_P_h
                                        * Referenced by: '<S214>/Proportional Gain'
                                        */
  real_T ProportionalAltitude_P_e;     /* Mask Parameter: ProportionalAltitude_P_e
                                        * Referenced by: '<S217>/Proportional Gain'
                                        */
  real_T PHeight_P_d;                  /* Mask Parameter: PHeight_P_d
                                        * Referenced by: '<S233>/Proportional Gain'
                                        */
  real_T PIDController1_P_a;           /* Mask Parameter: PIDController1_P_a
                                        * Referenced by: '<S235>/Proportional Gain'
                                        */
  real_T PIDController_P_p;            /* Mask Parameter: PIDController_P_p
                                        * Referenced by: '<S234>/Proportional Gain'
                                        */
  real_T PIDController3_P_a;           /* Mask Parameter: PIDController3_P_a
                                        * Referenced by: '<S237>/Proportional Gain'
                                        */
  real_T PIDController2_P_p;           /* Mask Parameter: PIDController2_P_p
                                        * Referenced by: '<S236>/Proportional Gain'
                                        */
  real_T PIDController5_P_b;           /* Mask Parameter: PIDController5_P_b
                                        * Referenced by: '<S239>/Proportional Gain'
                                        */
  real_T PIDController4_P_o;           /* Mask Parameter: PIDController4_P_o
                                        * Referenced by: '<S238>/Proportional Gain'
                                        */
  real_T AerodynamicForcesandMoments_S;/* Mask Parameter: AerodynamicForcesandMoments_S
                                        * Referenced by: '<S252>/reference area'
                                        */
  real_T DrydenWindTurbulenceModelDiscreteqr_T_on;/* Mask Parameter: DrydenWindTurbulenceModelDiscreteqr_T_on
                                                   * Referenced by:
                                                   *   '<S55>/Constant1'
                                                   *   '<S55>/Constant2'
                                                   *   '<S55>/Constant3'
                                                   *   '<S56>/Constant3'
                                                   */
  real_T DrydenWindTurbulenceModelDiscreteqr_TurbProb;/* Mask Parameter: DrydenWindTurbulenceModelDiscreteqr_TurbProb
                                                       * Referenced by: '<S73>/Probability of  Exceedance'
                                                       */
  real_T PIDController7_UpperSaturationLimit;/* Mask Parameter: PIDController7_UpperSaturationLimit
                                              * Referenced by: '<S176>/Saturate'
                                              */
  real_T PHeight_UpperSaturationLimit; /* Mask Parameter: PHeight_UpperSaturationLimit
                                        * Referenced by: '<S169>/Saturate'
                                        */
  real_T PIDController_UpperSaturationLimit;/* Mask Parameter: PIDController_UpperSaturationLimit
                                             * Referenced by: '<S170>/Saturate'
                                             */
  real_T ProportionalAltitude_UpperSaturationLimit;/* Mask Parameter: ProportionalAltitude_UpperSaturationLimit
                                                    * Referenced by: '<S189>/Saturate'
                                                    */
  real_T PHeight1_UpperSaturationLimit;/* Mask Parameter: PHeight1_UpperSaturationLimit
                                        * Referenced by: '<S201>/Saturate'
                                        */
  real_T PIDController4_UpperSaturationLimit;/* Mask Parameter: PIDController4_UpperSaturationLimit
                                              * Referenced by: '<S206>/Saturate'
                                              */
  real_T ProportionalAltitude_UpperSaturationLimit_o;/* Mask Parameter: ProportionalAltitude_UpperSaturationLimit_o
                                                      * Referenced by: '<S217>/Saturate'
                                                      */
  real_T PHeight_UpperSaturationLimit_g;/* Mask Parameter: PHeight_UpperSaturationLimit_g
                                         * Referenced by: '<S233>/Saturate'
                                         */
  real_T PIDController_UpperSaturationLimit_b;/* Mask Parameter: PIDController_UpperSaturationLimit_b
                                               * Referenced by: '<S234>/Saturate'
                                               */
  real_T CustomVariableMass6DOFEulerAnglespropio_Vm_0[3];/* Mask Parameter: CustomVariableMass6DOFEulerAnglespropio_Vm_0
                                                          * Referenced by: '<S243>/ub,vb,wb'
                                                          */
  real_T DrydenWindTurbulenceModelDiscreteqr_W20;/* Mask Parameter: DrydenWindTurbulenceModelDiscreteqr_W20
                                                  * Referenced by: '<S48>/Windspeed at 20ft (6m)'
                                                  */
  real_T WindShearModel_W_20;          /* Mask Parameter: WindShearModel_W_20
                                        * Referenced by: '<S49>/Wind speed at reference height'
                                        */
  real_T DrydenWindTurbulenceModelDiscreteqr_Wdeg;/* Mask Parameter: DrydenWindTurbulenceModelDiscreteqr_Wdeg
                                                   * Referenced by: '<S48>/Wind direction'
                                                   */
  real_T WindShearModel_Wdeg;          /* Mask Parameter: WindShearModel_Wdeg
                                        * Referenced by: '<S49>/Wind Direction'
                                        */
  real_T DrydenWindTurbulenceModelDiscreteqr_Wingspan;/* Mask Parameter: DrydenWindTurbulenceModelDiscreteqr_Wingspan
                                                       * Referenced by: '<S48>/Wingspan'
                                                       */
  real_T ThreeaxisInertialMeasurementUnit_a_bias[3];/* Mask Parameter: ThreeaxisInertialMeasurementUnit_a_bias
                                                     * Referenced by: '<S531>/Measurement bias'
                                                     */
  real_T ThreeaxisInertialMeasurementUnit_a_sf_cc[9];/* Mask Parameter: ThreeaxisInertialMeasurementUnit_a_sf_cc
                                                      * Referenced by: '<S531>/Scale factors & Cross-coupling  errors'
                                                      */
  real_T AerodynamicForcesandMoments_b;/* Mask Parameter: AerodynamicForcesandMoments_b
                                        * Referenced by: '<S252>/Constant'
                                        */
  real_T AerodynamicForcesandMoments_cbar;/* Mask Parameter: AerodynamicForcesandMoments_cbar
                                           * Referenced by: '<S252>/Constant1'
                                           */
  real_T CompareToConstant_const;      /* Mask Parameter: CompareToConstant_const
                                        * Referenced by: '<S325>/Constant'
                                        */
  real_T CompareToConstant_const_c;    /* Mask Parameter: CompareToConstant_const_c
                                        * Referenced by: '<S323>/Constant'
                                        */
  real_T CompareToConstant_const_o;    /* Mask Parameter: CompareToConstant_const_o
                                        * Referenced by: '<S326>/Constant'
                                        */
  real_T CompareToConstant_const_n;    /* Mask Parameter: CompareToConstant_const_n
                                        * Referenced by: '<S319>/Constant'
                                        */
  real_T CompareToConstant_const_p;    /* Mask Parameter: CompareToConstant_const_p
                                        * Referenced by: '<S317>/Constant'
                                        */
  real_T CompareToConstant_const_h;    /* Mask Parameter: CompareToConstant_const_h
                                        * Referenced by: '<S320>/Constant'
                                        */
  real_T CompareToConstant_const_i;    /* Mask Parameter: CompareToConstant_const_i
                                        * Referenced by: '<S517>/Constant'
                                        */
  real_T CompareToConstant_const_d;    /* Mask Parameter: CompareToConstant_const_d
                                        * Referenced by: '<S515>/Constant'
                                        */
  real_T CompareToConstant_const_n3;   /* Mask Parameter: CompareToConstant_const_n3
                                        * Referenced by: '<S518>/Constant'
                                        */
  real_T CompareToConstant_const_na;   /* Mask Parameter: CompareToConstant_const_na
                                        * Referenced by: '<S511>/Constant'
                                        */
  real_T CompareToConstant_const_j;    /* Mask Parameter: CompareToConstant_const_j
                                        * Referenced by: '<S509>/Constant'
                                        */
  real_T CompareToConstant_const_k;    /* Mask Parameter: CompareToConstant_const_k
                                        * Referenced by: '<S512>/Constant'
                                        */
  real_T Distanceintogustx_d_m;        /* Mask Parameter: Distanceintogustx_d_m
                                        * Referenced by: '<S50>/Distance into Gust (x) (Limited to gust length d)'
                                        */
  real_T Distanceintogusty_d_m;        /* Mask Parameter: Distanceintogusty_d_m
                                        * Referenced by: '<S47>/Distance into gust (y)'
                                        */
  real_T Distanceintogustz_d_m;        /* Mask Parameter: Distanceintogustz_d_m
                                        * Referenced by: '<S47>/Distance into gust (z)'
                                        */
  real_T DiscreteWindGustModel_d_m[3]; /* Mask Parameter: DiscreteWindGustModel_d_m
                                        * Referenced by: '<S47>/pi//Gust length'
                                        */
  real_T EstimateCenterofGravity_emass;/* Mask Parameter: EstimateCenterofGravity_emass
                                        * Referenced by: '<S336>/Constant'
                                        */
  real_T CustomVariableMass6DOFEulerAnglespropio_eul_0[3];/* Mask Parameter: CustomVariableMass6DOFEulerAnglespropio_eul_0
                                                           * Referenced by: '<S291>/phi theta psi'
                                                           */
  real_T EstimateCenterofGravity_fmass;/* Mask Parameter: EstimateCenterofGravity_fmass
                                        * Referenced by: '<S336>/Constant1'
                                        */
  real_T ThreeaxisInertialMeasurementUnit_g_bias[3];/* Mask Parameter: ThreeaxisInertialMeasurementUnit_g_bias
                                                     * Referenced by: '<S532>/Measurement bias'
                                                     */
  real_T ThreeaxisInertialMeasurementUnit_g_sens[3];/* Mask Parameter: ThreeaxisInertialMeasurementUnit_g_sens
                                                     * Referenced by: '<S532>/g-sensitive bias'
                                                     */
  real_T ThreeaxisInertialMeasurementUnit_g_sf_cc[9];/* Mask Parameter: ThreeaxisInertialMeasurementUnit_g_sf_cc
                                                      * Referenced by: '<S532>/Scale factors & Cross-coupling  errors '
                                                      */
  real_T ThreeaxisInertialMeasurementUnit_imu[3];/* Mask Parameter: ThreeaxisInertialMeasurementUnit_imu
                                                  * Referenced by: '<S531>/wl_ins'
                                                  */
  real_T JKFlipFlop_initial_condition; /* Mask Parameter: JKFlipFlop_initial_condition
                                        * Referenced by: '<S465>/J-K Flip-Flop'
                                        */
  real_T JKFlipFlop_initial_condition_g;/* Mask Parameter: JKFlipFlop_initial_condition_g
                                         * Referenced by: '<S466>/J-K Flip-Flop'
                                         */
  real_T JKFlipFlop_initial_condition_o;/* Mask Parameter: JKFlipFlop_initial_condition_o
                                         * Referenced by: '<S467>/J-K Flip-Flop'
                                         */
  real_T JKFlipFlop_initial_condition_c;/* Mask Parameter: JKFlipFlop_initial_condition_c
                                         * Referenced by: '<S468>/J-K Flip-Flop'
                                         */
  real_T JKFlipFlop_initial_condition_d;/* Mask Parameter: JKFlipFlop_initial_condition_d
                                         * Referenced by: '<S469>/J-K Flip-Flop'
                                         */
  real_T JKFlipFlop_initial_condition_k;/* Mask Parameter: JKFlipFlop_initial_condition_k
                                         * Referenced by: '<S470>/J-K Flip-Flop'
                                         */
  real_T JKFlipFlop_initial_condition_gr;/* Mask Parameter: JKFlipFlop_initial_condition_gr
                                          * Referenced by: '<S471>/J-K Flip-Flop'
                                          */
  real_T JKFlipFlop_initial_condition_gd;/* Mask Parameter: JKFlipFlop_initial_condition_gd
                                          * Referenced by: '<S472>/J-K Flip-Flop'
                                          */
  real_T JKFlipFlop_initial_condition_m;/* Mask Parameter: JKFlipFlop_initial_condition_m
                                         * Referenced by: '<S473>/J-K Flip-Flop'
                                         */
  real_T JKFlipFlop_initial_condition_l;/* Mask Parameter: JKFlipFlop_initial_condition_l
                                         * Referenced by: '<S474>/J-K Flip-Flop'
                                         */
  real_T JKFlipFlop_initial_condition_gh;/* Mask Parameter: JKFlipFlop_initial_condition_gh
                                          * Referenced by: '<S475>/J-K Flip-Flop'
                                          */
  real_T JKFlipFlop_initial_condition_de;/* Mask Parameter: JKFlipFlop_initial_condition_de
                                          * Referenced by: '<S476>/J-K Flip-Flop'
                                          */
  real_T InterpolateCG_matrix[6];      /* Mask Parameter: InterpolateCG_matrix
                                        * Referenced by:
                                        *   '<S340>/[0]'
                                        *   '<S340>/[1]'
                                        */
  real_T CheckAltitude_max;            /* Mask Parameter: CheckAltitude_max
                                        * Referenced by: '<S96>/max_val'
                                        */
  real_T CheckLatitude_max;            /* Mask Parameter: CheckLatitude_max
                                        * Referenced by: '<S97>/max_val'
                                        */
  real_T CheckLongitude_max;           /* Mask Parameter: CheckLongitude_max
                                        * Referenced by: '<S98>/max_val'
                                        */
  real_T Istimewithinmodellimits_max;  /* Mask Parameter: Istimewithinmodellimits_max
                                        * Referenced by: '<S100>/max_val'
                                        */
  real_T RotationMatrixtoVRMLRotation_maxzero;/* Mask Parameter: RotationMatrixtoVRMLRotation_maxzero
                                               * Referenced by:
                                               *   '<S558>/Switch'
                                               *   '<S559>/Dead Zone'
                                               */
  real_T CheckAltitude_min;            /* Mask Parameter: CheckAltitude_min
                                        * Referenced by: '<S96>/min_val'
                                        */
  real_T CheckLatitude_min;            /* Mask Parameter: CheckLatitude_min
                                        * Referenced by: '<S97>/min_val'
                                        */
  real_T CheckLongitude_min;           /* Mask Parameter: CheckLongitude_min
                                        * Referenced by: '<S98>/min_val'
                                        */
  real_T Istimewithinmodellimits_min;  /* Mask Parameter: Istimewithinmodellimits_min
                                        * Referenced by: '<S100>/min_val'
                                        */
  real_T CustomVariableMass6DOFEulerAnglespropio_pm_0[3];/* Mask Parameter: CustomVariableMass6DOFEulerAnglespropio_pm_0
                                                          * Referenced by: '<S243>/p,q,r '
                                                          */
  real_T FlatEarthtoLLA_psi;           /* Mask Parameter: FlatEarthtoLLA_psi
                                        * Referenced by: '<S313>/ref_pos'
                                        */
  real_T LLAtoFlatEarth_psi;           /* Mask Parameter: LLAtoFlatEarth_psi
                                        * Referenced by: '<S505>/ref_pos'
                                        */
  real_T DiscreteWindGustModel_t_0;    /* Mask Parameter: DiscreteWindGustModel_t_0
                                        * Referenced by: '<S47>/Gust start time'
                                        */
  real_T DiscreteWindGustModel_v_m[3]; /* Mask Parameter: DiscreteWindGustModel_v_m
                                        * Referenced by: '<S47>/Gust magnitude//2.0'
                                        */
  real_T CustomVariableMass6DOFEulerAnglespropio_xme_0[3];/* Mask Parameter: CustomVariableMass6DOFEulerAnglespropio_xme_0
                                                           * Referenced by: '<S243>/xe,ye,ze'
                                                           */
  int32_T CompareToConstant_const_f;   /* Mask Parameter: CompareToConstant_const_f
                                        * Referenced by: '<S228>/Constant'
                                        */
  int32_T CompareToConstant1_const;    /* Mask Parameter: CompareToConstant1_const
                                        * Referenced by: '<S229>/Constant'
                                        */
  int32_T CompareToConstant2_const;    /* Mask Parameter: CompareToConstant2_const
                                        * Referenced by: '<S230>/Constant'
                                        */
  int32_T CompareToConstant3_const;    /* Mask Parameter: CompareToConstant3_const
                                        * Referenced by: '<S231>/Constant'
                                        */
  int32_T CompareToConstant4_const;    /* Mask Parameter: CompareToConstant4_const
                                        * Referenced by: '<S232>/Constant'
                                        */
  ThrottleBus Throttle_Y0;             /* Computed Parameter: Throttle_Y0
                                        * Referenced by: '<S164>/Throttle'
                                        */
  ThrottleBus Throttle_Y0_b;           /* Computed Parameter: Throttle_Y0_b
                                        * Referenced by: '<S165>/Throttle'
                                        */
  ThrottleBus Throttle_Y0_e;           /* Computed Parameter: Throttle_Y0_e
                                        * Referenced by: '<S166>/Throttle'
                                        */
  ThrottleBus Throttle_Y0_m;           /* Computed Parameter: Throttle_Y0_m
                                        * Referenced by: '<S168>/Throttle'
                                        */
  ActuatorsBus Actuators_Y0;           /* Computed Parameter: Actuators_Y0
                                        * Referenced by: '<S163>/Actuators'
                                        */
  ActuatorsBus actuators_Y0;           /* Computed Parameter: actuators_Y0
                                        * Referenced by: '<S164>/actuators'
                                        */
  ActuatorsBus actuators_Y0_n;         /* Computed Parameter: actuators_Y0_n
                                        * Referenced by: '<S165>/actuators'
                                        */
  ActuatorsBus actuators_Y0_h;         /* Computed Parameter: actuators_Y0_h
                                        * Referenced by: '<S166>/actuators'
                                        */
  ActuatorsBus actuators_Y0_m;         /* Computed Parameter: actuators_Y0_m
                                        * Referenced by: '<S168>/actuators'
                                        */
  real_T uftinf_UpperSat;              /* Expression: inf
                                        * Referenced by: '<S49>/3ft-->inf'
                                        */
  real_T uftinf_LowerSat;              /* Expression: 3
                                        * Referenced by: '<S49>/3ft-->inf'
                                        */
  real_T hz0_Gain;                     /* Expression: 1/z0
                                        * Referenced by: '<S49>/h//z0'
                                        */
  real_T x_Y0;                         /* Expression: [0]
                                        * Referenced by: '<S50>/x'
                                        */
  real_T DistanceintoGustxLimitedtogustlengthd_IC;/* Expression: 0
                                                   * Referenced by: '<S50>/Distance into Gust (x) (Limited to gust length d)'
                                                   */
  real_T DistanceintoGustxLimitedtogustlengthd_LowerSat;/* Expression: 0
                                                         * Referenced by: '<S50>/Distance into Gust (x) (Limited to gust length d)'
                                                         */
  real_T pgw_Y0;                       /* Expression: 0
                                        * Referenced by: '<S67>/pgw'
                                        */
  real_T Constant2_Value;              /* Expression: 2.6
                                        * Referenced by: '<S67>/Constant2'
                                        */
  real_T u_Gain;                       /* Expression: 2*dt
                                        * Referenced by: '<S67>/2'
                                        */
  real_T Constant_Value;               /* Expression: 1
                                        * Referenced by: '<S67>/Constant'
                                        */
  real_T Constant1_Value;              /* Expression: 0.95
                                        * Referenced by: '<S67>/Constant1'
                                        */
  real_T Constant3_Value;              /* Expression: 1/3
                                        * Referenced by: '<S67>/Constant3'
                                        */
  real_T dt_Gain;                      /* Expression: dt
                                        * Referenced by: '<S67>/dt'
                                        */
  real_T UnitDelay_InitialCondition;   /* Expression: 0
                                        * Referenced by: '<S67>/Unit Delay'
                                        */
  real_T qgw_Y0;                       /* Expression: 0
                                        * Referenced by: '<S68>/qgw'
                                        */
  real_T Constant_Value_k;             /* Expression: 1
                                        * Referenced by: '<S68>/Constant'
                                        */
  real_T dt1_Gain;                     /* Expression: 4/pi
                                        * Referenced by: '<S68>/dt1'
                                        */
  real_T dt_Gain_j;                    /* Expression: dt
                                        * Referenced by: '<S68>/dt'
                                        */
  real_T UnitDelay_InitialCondition_f; /* Expression: 0
                                        * Referenced by: '<S68>/Unit Delay'
                                        */
  real_T UnitDelay1_InitialCondition;  /* Expression: 0
                                        * Referenced by: '<S68>/Unit Delay1'
                                        */
  real_T rgw_Y0;                       /* Expression: 0
                                        * Referenced by: '<S69>/rgw'
                                        */
  real_T Constant_Value_l;             /* Expression: 1
                                        * Referenced by: '<S69>/Constant'
                                        */
  real_T dt1_Gain_l;                   /* Expression: 3/pi
                                        * Referenced by: '<S69>/dt1'
                                        */
  real_T dt_Gain_k;                    /* Expression: dt
                                        * Referenced by: '<S69>/dt'
                                        */
  real_T UnitDelay_InitialCondition_k; /* Expression: 0
                                        * Referenced by: '<S69>/Unit Delay'
                                        */
  real_T UnitDelay1_InitialCondition_h;/* Expression: 0
                                        * Referenced by: '<S69>/Unit Delay1'
                                        */
  real_T ugw_Y0;                       /* Expression: 0
                                        * Referenced by: '<S70>/ugw'
                                        */
  real_T u_Gain_a;                     /* Expression: 2*dt
                                        * Referenced by: '<S70>/2'
                                        */
  real_T Constant_Value_b;             /* Expression: 1
                                        * Referenced by: '<S70>/Constant'
                                        */
  real_T dt_Gain_jf;                   /* Expression: dt
                                        * Referenced by: '<S70>/dt'
                                        */
  real_T UnitDelay_InitialCondition_m; /* Expression: 0
                                        * Referenced by: '<S70>/Unit Delay'
                                        */
  real_T vgw_Y0;                       /* Expression: 0
                                        * Referenced by: '<S71>/vgw'
                                        */
  real_T u_Gain_g;                     /* Expression: 2*dt
                                        * Referenced by: '<S71>/2'
                                        */
  real_T Constant_Value_p;             /* Expression: 1
                                        * Referenced by: '<S71>/Constant'
                                        */
  real_T dt_Gain_f;                    /* Expression: dt
                                        * Referenced by: '<S71>/dt'
                                        */
  real_T UnitDelay_InitialCondition_e; /* Expression: 0
                                        * Referenced by: '<S71>/Unit Delay'
                                        */
  real_T wgw_Y0;                       /* Expression: 0
                                        * Referenced by: '<S72>/wgw'
                                        */
  real_T u_Gain_h;                     /* Expression: 2*dt
                                        * Referenced by: '<S72>/2'
                                        */
  real_T Constant_Value_o;             /* Expression: 1
                                        * Referenced by: '<S72>/Constant'
                                        */
  real_T dt_Gain_o;                    /* Expression: dt
                                        * Referenced by: '<S72>/dt'
                                        */
  real_T UnitDelay_InitialCondition_kv;/* Expression: 0
                                        * Referenced by: '<S72>/Unit Delay'
                                        */
  real_T Gain_Gain;                    /* Expression: 1
                                        * Referenced by: '<S77>/Gain'
                                        */
  real_T max_height_low_Value;         /* Expression: max_height_low
                                        * Referenced by: '<S75>/max_height_low'
                                        */
  real_T min_height_high_Value;        /* Expression: min_height_high
                                        * Referenced by: '<S75>/min_height_high'
                                        */
  real_T Gain_Gain_c;                  /* Expression: 1
                                        * Referenced by: '<S85>/Gain'
                                        */
  real_T max_height_low_Value_i;       /* Expression: max_height_low
                                        * Referenced by: '<S83>/max_height_low'
                                        */
  real_T min_height_high_Value_n;      /* Expression: min_height_high
                                        * Referenced by: '<S83>/min_height_high'
                                        */
  real_T pp13_Y0[13];                  /* Expression: ones(1,maxdef+1)
                                        * Referenced by: '<S122>/pp[13]'
                                        */
  real_T Constant_Value_k0;            /* Expression: 1
                                        * Referenced by: '<S122>/Constant'
                                        */
  real_T pp13_Y0_d[13];                /* Expression: ones(1,maxdef+1)
                                        * Referenced by: '<S123>/pp[13]'
                                        */
  real_T k1313_Value[169];             /* Expression: k
                                        * Referenced by: '<S123>/k[13][13]'
                                        */
  real_T bpp_Y0;                       /* Expression: 0
                                        * Referenced by: '<S118>/bpp'
                                        */
  real_T Constant_Value_g;             /* Expression: fm(2)
                                        * Referenced by: '<S118>/Constant'
                                        */
  real_T Constant1_Value_f;            /* Expression: 1
                                        * Referenced by: '<S118>/Constant1'
                                        */
  real_T UnitDelay1_InitialCondition_l[13];/* Expression: ones(1,maxdef+1)
                                            * Referenced by: '<S118>/Unit Delay1'
                                            */
  real_T Constant_Value_lr;            /* Expression: 1
                                        * Referenced by: '<S126>/Constant'
                                        */
  real_T Gain1_Gain;                   /* Expression: 1
                                        * Referenced by: '<S126>/Gain1'
                                        */
  real_T Gain2_Gain;                   /* Expression: 1
                                        * Referenced by: '<S126>/Gain2'
                                        */
  real_T Constant_Value_n;             /* Expression: 1
                                        * Referenced by: '<S128>/Constant'
                                        */
  real_T Constant_Value_km;            /* Expression: 1
                                        * Referenced by: '<S129>/Constant'
                                        */
  real_T Constant1_Value_d;            /* Expression: 0
                                        * Referenced by: '<S132>/Constant1'
                                        */
  real_T Constant_Value_i;             /* Expression: 0
                                        * Referenced by: '<S132>/Constant'
                                        */
  real_T Switch_Threshold;             /* Expression: 0.5
                                        * Referenced by: '<S132>/Switch'
                                        */
  real_T k1313_Value_m[169];           /* Expression: k
                                        * Referenced by: '<S132>/k[13][13]'
                                        */
  real_T Switch1_Threshold;            /* Expression: 0.5
                                        * Referenced by: '<S132>/Switch1'
                                        */
  real_T dp1313_Y0[169];               /* Expression: zeros(maxdef+1,maxdef+1)
                                        * Referenced by: '<S116>/dp[13][13]'
                                        */
  real_T snorm169_Y0[169];             /* Expression: snorm
                                        * Referenced by: '<S116>/snorm[169]'
                                        */
  real_T UnitDelay_InitialCondition_h[169];/* Expression: zeros(maxdef+1,maxdef+1)
                                            * Referenced by: '<S116>/Unit Delay'
                                            */
  real_T UnitDelay1_InitialCondition_g[169];/* Expression: snorm
                                             * Referenced by: '<S116>/Unit Delay1'
                                             */
  real_T Merge1_InitialOutput;         /* Expression: 0
                                        * Referenced by: '<S116>/Merge1'
                                        */
  real_T Merge_InitialOutput;          /* Expression: 0
                                        * Referenced by: '<S116>/Merge'
                                        */
  real_T zerosmaxdef1maxdef1_Value[169];/* Expression: zeros(maxdef+1,maxdef+1)
                                         * Referenced by: '<S144>/zeros(maxdef+1,maxdef+1)'
                                         */
  real_T Gain_Gain_d;                  /* Expression: 1
                                        * Referenced by: '<S145>/Gain'
                                        */
  real_T tc1313_Y0[169];               /* Expression: zeros(maxdef+1,maxdef+1)
                                        * Referenced by: '<S117>/tc[13][13]'
                                        */
  real_T UnitDelay_InitialCondition_e2[169];/* Expression: zeros(maxdef+1,maxdef+1)
                                             * Referenced by: '<S117>/Unit Delay'
                                             */
  real_T cmaxdefmaxdef_Value[169];     /* Expression: c
                                        * Referenced by: '<S117>/c[maxdef][maxdef]'
                                        */
  real_T cdmaxdefmaxdef_Value[169];    /* Expression: cd
                                        * Referenced by: '<S117>/cd[maxdef][maxdef]'
                                        */
  real_T UnitDelay_InitialCondition_o[169];/* Expression: zeros(maxdef+1,maxdef+1)
                                            * Referenced by: '<S144>/Unit Delay'
                                            */
  real_T bt_Y0;                        /* Expression: 0
                                        * Referenced by: '<S114>/bt'
                                        */
  real_T bp_Y0;                        /* Expression: 0
                                        * Referenced by: '<S114>/bp'
                                        */
  real_T br_Y0;                        /* Expression: 0
                                        * Referenced by: '<S114>/br'
                                        */
  real_T bpp_Y0_f;                     /* Expression: 0
                                        * Referenced by: '<S114>/bpp'
                                        */
  real_T Constant1_Value_p;            /* Expression: 1
                                        * Referenced by: '<S120>/Constant1'
                                        */
  real_T Merge_InitialOutput_d;        /* Expression: 0
                                        * Referenced by: '<S120>/Merge'
                                        */
  real_T fm_Value[13];                 /* Expression: fm
                                        * Referenced by: '<S115>/fm'
                                        */
  real_T Merge1_InitialOutput_a;       /* Expression: 0
                                        * Referenced by: '<S120>/Merge1'
                                        */
  real_T fn_Value[13];                 /* Expression: fn
                                        * Referenced by: '<S115>/fn'
                                        */
  real_T Constant1_Value_o;            /* Expression: 0
                                        * Referenced by: '<S121>/Constant1'
                                        */
  real_T UnitDelay1_InitialCondition_a;/* Expression: 0
                                        * Referenced by: '<S115>/Unit Delay1'
                                        */
  real_T UnitDelay3_InitialCondition;  /* Expression: 0
                                        * Referenced by: '<S115>/Unit Delay3'
                                        */
  real_T UnitDelay2_InitialCondition;  /* Expression: 0
                                        * Referenced by: '<S115>/Unit Delay2'
                                        */
  real_T UnitDelay4_InitialCondition;  /* Expression: 0
                                        * Referenced by: '<S115>/Unit Delay4'
                                        */
  real_T btbpbrbpp_Y0[4];              /* Expression: [0 0 0 0]
                                        * Referenced by: '<S106>/bt,bp,br,bpp'
                                        */
  real_T UnitDelay_InitialCondition_c; /* Expression: 0
                                        * Referenced by: '<S106>/Unit Delay'
                                        */
  real_T UnitDelay2_InitialCondition_j[4];/* Expression: [0 0 0 0]
                                           * Referenced by: '<S106>/Unit Delay2'
                                           */
  real_T r_Y0;                         /* Expression: 6378.137
                                        * Referenced by: '<S107>/r'
                                        */
  real_T ct_Y0;                        /* Expression: 1
                                        * Referenced by: '<S107>/ct'
                                        */
  real_T st_Y0;                        /* Expression: 0
                                        * Referenced by: '<S107>/st'
                                        */
  real_T sa_Y0;                        /* Expression: 0
                                        * Referenced by: '<S107>/sa'
                                        */
  real_T ca_Y0;                        /* Expression: 0
                                        * Referenced by: '<S107>/ca'
                                        */
  real_T b_Value;                      /* Expression: 6356.7523142
                                        * Referenced by: '<S107>/b'
                                        */
  real_T a_Value;                      /* Expression: 6378.137
                                        * Referenced by: '<S107>/a'
                                        */
  real_T Gain_Gain_f;                  /* Expression: 2
                                        * Referenced by: '<S152>/Gain'
                                        */
  real_T Constant_Value_j;             /* Expression: 1
                                        * Referenced by: '<S154>/Constant'
                                        */
  real_T sp11_Y0[11];                  /* Expression: (1:(maxdef-1))
                                        * Referenced by: '<S155>/sp[11]'
                                        */
  real_T cp11_Y0[11];                  /* Expression: (1:(maxdef-1))
                                        * Referenced by: '<S155>/cp[11]'
                                        */
  real_T Constant_Value_m[11];         /* Expression: [1:maxdef-1]
                                        * Referenced by: '<S155>/Constant'
                                        */
  real_T UnitDelay1_InitialCondition_h1;/* Expression: 0
                                         * Referenced by: '<S155>/Unit Delay1'
                                         */
  real_T ForIterator_IterationLimit;   /* Expression: maxdef-1
                                        * Referenced by: '<S155>/For Iterator'
                                        */
  real_T cpm1spm1_Threshold;           /* Expression: 1
                                        * Referenced by: '<S155>/cp[m-1] sp[m-1]'
                                        */
  real_T Constant1_Value_l[11];        /* Expression: [1:maxdef-1]
                                        * Referenced by: '<S155>/Constant1'
                                        */
  real_T sp13_Y0[13];                  /* Expression: [0 0 (1:(maxdef-1))]
                                        * Referenced by: '<S108>/sp[13]'
                                        */
  real_T cp13_Y0[13];                  /* Expression: [1 1 (1:(maxdef-1))]
                                        * Referenced by: '<S108>/cp[13]'
                                        */
  real_T Gain_Gain_o;                  /* Expression: 1
                                        * Referenced by: '<S108>/Gain'
                                        */
  real_T Gain1_Gain_h;                 /* Expression: 1
                                        * Referenced by: '<S108>/Gain1'
                                        */
  real_T cp1_Value;                    /* Expression: 1
                                        * Referenced by: '<S108>/cp[1]'
                                        */
  real_T sp1_Value;                    /* Expression: 0
                                        * Referenced by: '<S108>/sp[1]'
                                        */
  real_T Constant9_Value;              /* Expression: 30
                                        * Referenced by: '<S163>/Constant9'
                                        */
  real_T Gain_Gain_i;                  /* Expression: pi/180
                                        * Referenced by: '<S163>/Gain'
                                        */
  real_T Constant10_Value;             /* Expression: 30
                                        * Referenced by: '<S163>/Constant10'
                                        */
  real_T Gain1_Gain_i;                 /* Expression: pi/180
                                        * Referenced by: '<S163>/Gain1'
                                        */
  real_T Constant1_Value_b;            /* Expression: 30
                                        * Referenced by: '<S163>/Constant1'
                                        */
  real_T Gain2_Gain_g;                 /* Expression: pi/180
                                        * Referenced by: '<S163>/Gain2'
                                        */
  real_T Constant2_Value_h;            /* Expression: 30
                                        * Referenced by: '<S163>/Constant2'
                                        */
  real_T Gain3_Gain;                   /* Expression: pi/180
                                        * Referenced by: '<S163>/Gain3'
                                        */
  real_T Constant21_Value;             /* Expression: 1
                                        * Referenced by: '<S163>/Constant21'
                                        */
  real_T IC_Value;                     /* Expression: 0
                                        * Referenced by: '<S163>/IC'
                                        */
  real_T Constant22_Value;             /* Expression: 1
                                        * Referenced by: '<S163>/Constant22'
                                        */
  real_T IC1_Value;                    /* Expression: 0
                                        * Referenced by: '<S163>/IC1'
                                        */
  real_T Constant23_Value;             /* Expression: 1
                                        * Referenced by: '<S163>/Constant23'
                                        */
  real_T IC2_Value;                    /* Expression: 0
                                        * Referenced by: '<S163>/IC2'
                                        */
  real_T Constant24_Value;             /* Expression: 1
                                        * Referenced by: '<S163>/Constant24'
                                        */
  real_T IC3_Value;                    /* Expression: 0
                                        * Referenced by: '<S163>/IC3'
                                        */
  real_T Constant20_Value;             /* Expression: 1
                                        * Referenced by: '<S163>/Constant20'
                                        */
  real_T IC4_Value;                    /* Expression: 0
                                        * Referenced by: '<S163>/IC4'
                                        */
  real_T Gain8_Gain;                   /* Expression: -1
                                        * Referenced by: '<S164>/Gain8'
                                        */
  real_T Gain_Gain_b;                  /* Expression: pi/180
                                        * Referenced by: '<S164>/Gain'
                                        */
  real_T Gain1_Gain_p;                 /* Expression: pi/180
                                        * Referenced by: '<S164>/Gain1'
                                        */
  real_T Gain2_Gain_m;                 /* Expression: pi/180
                                        * Referenced by: '<S164>/Gain2'
                                        */
  real_T Gain3_Gain_d;                 /* Expression: pi/180
                                        * Referenced by: '<S164>/Gain3'
                                        */
  real_T Constant18_Value;             /* Expression: 19.5
                                        * Referenced by: '<S164>/Constant18'
                                        */
  real_T Integrator_IC;                /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S176>/Integrator'
                                        */
  real_T Filter_IC;                    /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S176>/Filter'
                                        */
  real_T Constant15_Value;             /* Expression: 20
                                        * Referenced by: '<S164>/Constant15'
                                        */
  real_T Saturation_UpperSat;          /* Expression: 25
                                        * Referenced by: '<S164>/Saturation'
                                        */
  real_T Saturation_LowerSat;          /* Expression: -25
                                        * Referenced by: '<S164>/Saturation'
                                        */
  real_T Integrator_IC_o;              /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S169>/Integrator'
                                        */
  real_T Filter_IC_m;                  /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S169>/Filter'
                                        */
  real_T Integrator_IC_d;              /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S171>/Integrator'
                                        */
  real_T Filter_IC_a;                  /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S171>/Filter'
                                        */
  real_T Integrator_IC_j;              /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S170>/Integrator'
                                        */
  real_T Filter_IC_o;                  /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S170>/Filter'
                                        */
  real_T Gain4_Gain;                   /* Expression: pi/180
                                        * Referenced by: '<S164>/Gain4'
                                        */
  real_T Constant14_Value;             /* Expression: 0
                                        * Referenced by: '<S164>/Constant14'
                                        */
  real_T Integrator_IC_e;              /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S173>/Integrator'
                                        */
  real_T Filter_IC_ar;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S173>/Filter'
                                        */
  real_T Integrator_IC_n;              /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S172>/Integrator'
                                        */
  real_T Filter_IC_f;                  /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S172>/Filter'
                                        */
  real_T Constant11_Value;             /* Expression: 0
                                        * Referenced by: '<S164>/Constant11'
                                        */
  real_T Gain5_Gain;                   /* Expression: pi/180
                                        * Referenced by: '<S164>/Gain5'
                                        */
  real_T Constant16_Value;             /* Expression: 0
                                        * Referenced by: '<S164>/Constant16'
                                        */
  real_T Integrator_IC_f;              /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S175>/Integrator'
                                        */
  real_T Filter_IC_e;                  /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S175>/Filter'
                                        */
  real_T Integrator_IC_fq;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S174>/Integrator'
                                        */
  real_T Filter_IC_c;                  /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S174>/Filter'
                                        */
  real_T Constant17_Value;             /* Expression: 10*pi/180
                                        * Referenced by: '<S164>/Constant17'
                                        */
  real_T Constant12_Value;             /* Expression: 20*pi/180
                                        * Referenced by: '<S164>/Constant12'
                                        */
  real_T Gain6_Gain;                   /* Expression: pi/180
                                        * Referenced by: '<S164>/Gain6'
                                        */
  real_T Constant13_Value;             /* Expression: 20*pi/180
                                        * Referenced by: '<S164>/Constant13'
                                        */
  real_T Gain7_Gain;                   /* Expression: pi/180
                                        * Referenced by: '<S164>/Gain7'
                                        */
  real_T Constant6_Value[2];           /* Expression: [0,0]
                                        * Referenced by: '<S164>/Constant6'
                                        */
  real_T Constant4_Value;              /* Expression: 1
                                        * Referenced by: '<S164>/Constant4'
                                        */
  real_T Constant5_Value;              /* Expression: 0
                                        * Referenced by: '<S164>/Constant5'
                                        */
  real_T Constant7_Value;              /* Expression: 0
                                        * Referenced by: '<S164>/Constant7'
                                        */
  real_T Constant8_Value;              /* Expression: 0
                                        * Referenced by: '<S164>/Constant8'
                                        */
  real_T Constant3_Value_m;            /* Expression: 0
                                        * Referenced by: '<S164>/Constant3'
                                        */
  real_T Gain_Gain_dd;                 /* Expression: pi/180
                                        * Referenced by: '<S165>/Gain'
                                        */
  real_T Gain1_Gain_e;                 /* Expression: pi/180
                                        * Referenced by: '<S165>/Gain1'
                                        */
  real_T Gain2_Gain_mx;                /* Expression: pi/180
                                        * Referenced by: '<S165>/Gain2'
                                        */
  real_T Gain3_Gain_e;                 /* Expression: pi/180
                                        * Referenced by: '<S165>/Gain3'
                                        */
  real_T Constant21_Value_f;           /* Expression: 0
                                        * Referenced by: '<S165>/Constant21'
                                        */
  real_T Constant22_Value_i;           /* Expression: 0
                                        * Referenced by: '<S165>/Constant22'
                                        */
  real_T Constant23_Value_p;           /* Expression: 0
                                        * Referenced by: '<S165>/Constant23'
                                        */
  real_T Constant24_Value_j;           /* Expression: 0
                                        * Referenced by: '<S165>/Constant24'
                                        */
  real_T Constant20_Value_j;           /* Expression: 0
                                        * Referenced by: '<S165>/Constant20'
                                        */
  real_T Constant_Value_d;             /* Expression: 13
                                        * Referenced by: '<S165>/Constant'
                                        */
  real_T Constant13_Value_f;           /* Expression: 0
                                        * Referenced by: '<S165>/Constant13'
                                        */
  real_T Constant4_Value_d;            /* Expression: 45*pi/180
                                        * Referenced by: '<S165>/Constant4'
                                        */
  real_T Constant6_Value_m[2];         /* Expression: [0,0]
                                        * Referenced by: '<S165>/Constant6'
                                        */
  real_T Filter_IC_j;                  /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S179>/Filter'
                                        */
  real_T Integrator_IC_i;              /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S179>/Integrator'
                                        */
  real_T LateralPositionCmd_Value;     /* Expression: 0
                                        * Referenced by: '<S165>/Lateral Position Cmd'
                                        */
  real_T Saturation9_UpperSat;         /* Expression: 10
                                        * Referenced by: '<S165>/Saturation9'
                                        */
  real_T Saturation9_LowerSat;         /* Expression: -10
                                        * Referenced by: '<S165>/Saturation9'
                                        */
  real_T Integrator_IC_o0;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S184>/Integrator'
                                        */
  real_T Filter_IC_eq;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S184>/Filter'
                                        */
  real_T Saturation3_UpperSat;         /* Expression: 1.5*9.81
                                        * Referenced by: '<S165>/Saturation3'
                                        */
  real_T Saturation3_LowerSat;         /* Expression: 0.5*9.81
                                        * Referenced by: '<S165>/Saturation3'
                                        */
  real_T Integrator_IC_p;              /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S180>/Integrator'
                                        */
  real_T Filter_IC_b;                  /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S180>/Filter'
                                        */
  real_T Constant13_Value_j;           /* Expression: 0
                                        * Referenced by: '<S193>/Constant13'
                                        */
  real_T Integrator_IC_jo;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S185>/Integrator'
                                        */
  real_T Filter_IC_eg;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S185>/Filter'
                                        */
  real_T Saturation2_UpperSat;         /* Expression: 5
                                        * Referenced by: '<S165>/Saturation2'
                                        */
  real_T Saturation2_LowerSat;         /* Expression: -5
                                        * Referenced by: '<S165>/Saturation2'
                                        */
  real_T Integrator_IC_k;              /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S183>/Integrator'
                                        */
  real_T Filter_IC_h;                  /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S183>/Filter'
                                        */
  real_T Saturation1_UpperSat;         /* Expression: 25*pi/180
                                        * Referenced by: '<S165>/Saturation1'
                                        */
  real_T Saturation1_LowerSat;         /* Expression: -25*pi/180
                                        * Referenced by: '<S165>/Saturation1'
                                        */
  real_T Integrator_IC_jb;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S191>/Integrator'
                                        */
  real_T Filter_IC_f1;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S191>/Filter'
                                        */
  real_T Integrator_IC_m;              /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S187>/Integrator'
                                        */
  real_T Filter_IC_cm;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S187>/Filter'
                                        */
  real_T dRoll_Gain;                   /* Expression: LD.FCS.Quadcopter.Roll.dRoll
                                        * Referenced by: '<S197>/dRoll'
                                        */
  real_T Gain_Gain_di;                 /* Expression: -1
                                        * Referenced by: '<S197>/Gain'
                                        */
  real_T dThrottle_Gain;               /* Expression: LD.FCS.Quadcopter.Throttle.dThrottle
                                        * Referenced by: '<S195>/dThrottle'
                                        */
  real_T Bias_Bias;                    /* Expression: LD.FCS.Quadcopter.Throttle.bias
                                        * Referenced by: '<S195>/Bias'
                                        */
  real_T Integrator_IC_jt;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S192>/Integrator'
                                        */
  real_T Filter_IC_hb;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S192>/Filter'
                                        */
  real_T Integrator_IC_p2;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S188>/Integrator'
                                        */
  real_T Filter_IC_or;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S188>/Filter'
                                        */
  real_T dYaw_Gain;                    /* Expression: LD.FCS.Quadcopter.Yaw.dYaw
                                        * Referenced by: '<S198>/dYaw'
                                        */
  real_T Saturation7_UpperSat;         /* Expression: 10
                                        * Referenced by: '<S165>/Saturation7'
                                        */
  real_T Saturation7_LowerSat;         /* Expression: -10
                                        * Referenced by: '<S165>/Saturation7'
                                        */
  real_T Integrator_IC_ir;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S182>/Integrator'
                                        */
  real_T Filter_IC_p;                  /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S182>/Filter'
                                        */
  real_T Saturation6_UpperSat;         /* Expression: 25*pi/180
                                        * Referenced by: '<S165>/Saturation6'
                                        */
  real_T Saturation6_LowerSat;         /* Expression: -25*pi/180
                                        * Referenced by: '<S165>/Saturation6'
                                        */
  real_T Integrator_IC_fd;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S190>/Integrator'
                                        */
  real_T Filter_IC_ao;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S190>/Filter'
                                        */
  real_T Integrator_IC_ev;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S186>/Integrator'
                                        */
  real_T Filter_IC_i;                  /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S186>/Filter'
                                        */
  real_T dPitch_Gain;                  /* Expression: LD.FCS.Quadcopter.Pitch.dPitch
                                        * Referenced by: '<S196>/dPitch'
                                        */
  real_T Saturation_UpperSat_p;        /* Expression: 1
                                        * Referenced by: '<S193>/Saturation'
                                        */
  real_T Saturation_LowerSat_c;        /* Expression: 0
                                        * Referenced by: '<S193>/Saturation'
                                        */
  real_T motor2position_Value[3];      /* Expression: LD.Propulsion.motor2.Position
                                        * Referenced by: '<S196>/motor2 position '
                                        */
  real_T Constant1_Value_n[3];         /* Expression: LD.Inertia.CG
                                        * Referenced by: '<S196>/Constant1'
                                        */
  real_T motor2position1_Value[3];     /* Expression: LD.Propulsion.motor3.Position
                                        * Referenced by: '<S196>/motor2 position 1'
                                        */
  real_T Constant2_Value_d[3];         /* Expression: LD.Inertia.CG
                                        * Referenced by: '<S196>/Constant2'
                                        */
  real_T Gain_Gain_g;                  /* Expression: -1
                                        * Referenced by: '<S198>/Gain'
                                        */
  real_T Gain1_Gain_e4;                /* Expression: -1
                                        * Referenced by: '<S197>/Gain1'
                                        */
  real_T Saturation1_UpperSat_h;       /* Expression: 1
                                        * Referenced by: '<S193>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_i;       /* Expression: 0
                                        * Referenced by: '<S193>/Saturation1'
                                        */
  real_T Saturation2_UpperSat_m;       /* Expression: 1
                                        * Referenced by: '<S193>/Saturation2'
                                        */
  real_T Saturation2_LowerSat_k;       /* Expression: 0
                                        * Referenced by: '<S193>/Saturation2'
                                        */
  real_T Gain1_Gain_k;                 /* Expression: -1
                                        * Referenced by: '<S198>/Gain1'
                                        */
  real_T Saturation3_UpperSat_p;       /* Expression: 1
                                        * Referenced by: '<S193>/Saturation3'
                                        */
  real_T Saturation3_LowerSat_c;       /* Expression: 0
                                        * Referenced by: '<S193>/Saturation3'
                                        */
  real_T Integrator_IC_ph;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S189>/Integrator'
                                        */
  real_T Filter_IC_pe;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S189>/Filter'
                                        */
  real_T Filter_IC_g;                  /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S181>/Filter'
                                        */
  real_T Integrator_IC_nq;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S181>/Integrator'
                                        */
  real_T Saturation6_UpperSat_i;       /* Expression: 25*pi/180
                                        * Referenced by: '<S166>/Saturation6'
                                        */
  real_T Saturation6_LowerSat_f;       /* Expression: -25*pi/180
                                        * Referenced by: '<S166>/Saturation6'
                                        */
  real_T Gain8_Gain_d;                 /* Expression: -1
                                        * Referenced by: '<S166>/Gain8'
                                        */
  real_T Constant2_Value_o;            /* Expression: 20*pi/180
                                        * Referenced by: '<S166>/Constant2'
                                        */
  real_T Constant25_Value;             /* Expression: 0
                                        * Referenced by: '<S166>/Constant25'
                                        */
  real_T Constant21_Value_fz;          /* Expression: 1
                                        * Referenced by: '<S166>/Constant21'
                                        */
  real_T Constant22_Value_f;           /* Expression: 0
                                        * Referenced by: '<S166>/Constant22'
                                        */
  real_T Constant23_Value_f;           /* Expression: 0
                                        * Referenced by: '<S166>/Constant23'
                                        */
  real_T Constant24_Value_e;           /* Expression: 0
                                        * Referenced by: '<S166>/Constant24'
                                        */
  real_T Constant20_Value_e;           /* Expression: 0
                                        * Referenced by: '<S166>/Constant20'
                                        */
  real_T Constant18_Value_e;           /* Expression: 100
                                        * Referenced by: '<S166>/Constant18'
                                        */
  real_T Constant6_Value_j[2];         /* Expression: [0,0]
                                        * Referenced by: '<S166>/Constant6'
                                        */
  real_T Saturation11_UpperSat;        /* Expression: 25
                                        * Referenced by: '<S166>/Saturation11'
                                        */
  real_T Saturation11_LowerSat;        /* Expression: -25
                                        * Referenced by: '<S166>/Saturation11'
                                        */
  real_T Integrator_IC_kg;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S201>/Integrator'
                                        */
  real_T Filter_IC_l;                  /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S201>/Filter'
                                        */
  real_T Integrator_IC_ig;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S207>/Integrator'
                                        */
  real_T Filter_IC_pt;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S207>/Filter'
                                        */
  real_T Integrator_IC_nu;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S206>/Integrator'
                                        */
  real_T Filter_IC_ms;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S206>/Filter'
                                        */
  real_T Gain4_Gain_f;                 /* Expression: pi/180
                                        * Referenced by: '<S166>/Gain4'
                                        */
  real_T Constant17_Value_n;           /* Expression: 0
                                        * Referenced by: '<S166>/Constant17'
                                        */
  real_T Integrator_IC_l;              /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S209>/Integrator'
                                        */
  real_T Filter_IC_em;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S209>/Filter'
                                        */
  real_T Integrator_IC_np;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S208>/Integrator'
                                        */
  real_T Filter_IC_gj;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S208>/Filter'
                                        */
  real_T Constant8_Value_l;            /* Expression: 0
                                        * Referenced by: '<S166>/Constant8'
                                        */
  real_T Gain5_Gain_p;                 /* Expression: pi/180
                                        * Referenced by: '<S166>/Gain5'
                                        */
  real_T Constant19_Value;             /* Expression: 0
                                        * Referenced by: '<S166>/Constant19'
                                        */
  real_T Integrator_IC_ot;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S211>/Integrator'
                                        */
  real_T Filter_IC_lx;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S211>/Filter'
                                        */
  real_T Integrator_IC_b;              /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S210>/Integrator'
                                        */
  real_T Filter_IC_gm;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S210>/Filter'
                                        */
  real_T Constant14_Value_p;           /* Expression: 20*pi/180
                                        * Referenced by: '<S166>/Constant14'
                                        */
  real_T Gain6_Gain_h;                 /* Expression: pi/180
                                        * Referenced by: '<S166>/Gain6'
                                        */
  real_T Constant16_Value_n;           /* Expression: 20*pi/180
                                        * Referenced by: '<S166>/Constant16'
                                        */
  real_T Gain7_Gain_f;                 /* Expression: pi/180
                                        * Referenced by: '<S166>/Gain7'
                                        */
  real_T Constant_Value_lj;            /* Expression: 13
                                        * Referenced by: '<S166>/Constant'
                                        */
  real_T Constant1_Value_pn;           /* Expression: 0*pi/180
                                        * Referenced by: '<S166>/Constant1'
                                        */
  real_T Constant13_Value_m;           /* Expression: 0
                                        * Referenced by: '<S166>/Constant13'
                                        */
  real_T Constant3_Value_n;            /* Expression: 1
                                        * Referenced by: '<S166>/Constant3'
                                        */
  real_T Constant4_Value_b;            /* Expression: 45*pi/180
                                        * Referenced by: '<S166>/Constant4'
                                        */
  real_T Filter_IC_is;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S200>/Filter'
                                        */
  real_T Integrator_IC_br;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S200>/Integrator'
                                        */
  real_T LateralPositionCmd_Value_c;   /* Expression: 0
                                        * Referenced by: '<S166>/Lateral Position Cmd'
                                        */
  real_T Saturation9_UpperSat_i;       /* Expression: 10
                                        * Referenced by: '<S166>/Saturation9'
                                        */
  real_T Saturation9_LowerSat_j;       /* Expression: -10
                                        * Referenced by: '<S166>/Saturation9'
                                        */
  real_T Integrator_IC_n5;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S212>/Integrator'
                                        */
  real_T Filter_IC_gn;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S212>/Filter'
                                        */
  real_T Saturation3_UpperSat_n;       /* Expression: 1.5*9.81
                                        * Referenced by: '<S166>/Saturation3'
                                        */
  real_T Saturation3_LowerSat_e;       /* Expression: 0.5*9.81
                                        * Referenced by: '<S166>/Saturation3'
                                        */
  real_T Integrator_IC_fx;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S202>/Integrator'
                                        */
  real_T Filter_IC_d;                  /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S202>/Filter'
                                        */
  real_T Integrator_IC_a;              /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S213>/Integrator'
                                        */
  real_T Filter_IC_mu;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S213>/Filter'
                                        */
  real_T Saturation2_UpperSat_d;       /* Expression: 5
                                        * Referenced by: '<S166>/Saturation2'
                                        */
  real_T Saturation2_LowerSat_d;       /* Expression: -5
                                        * Referenced by: '<S166>/Saturation2'
                                        */
  real_T Integrator_IC_ol;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S205>/Integrator'
                                        */
  real_T Filter_IC_jf;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S205>/Filter'
                                        */
  real_T Saturation1_UpperSat_o;       /* Expression: 25*pi/180
                                        * Referenced by: '<S166>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_p;       /* Expression: -25*pi/180
                                        * Referenced by: '<S166>/Saturation1'
                                        */
  real_T Integrator_IC_c;              /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S219>/Integrator'
                                        */
  real_T Filter_IC_og;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S219>/Filter'
                                        */
  real_T Integrator_IC_l5;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S215>/Integrator'
                                        */
  real_T Filter_IC_gmg;                /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S215>/Filter'
                                        */
  real_T dRoll_Gain_g;                 /* Expression: LD.FCS.Quadcopter.Roll.dRoll
                                        * Referenced by: '<S225>/dRoll'
                                        */
  real_T Gain_Gain_ci;                 /* Expression: -1
                                        * Referenced by: '<S225>/Gain'
                                        */
  real_T dThrottle_Gain_l;             /* Expression: LD.FCS.Quadcopter.Throttle.dThrottle
                                        * Referenced by: '<S223>/dThrottle'
                                        */
  real_T Bias_Bias_o;                  /* Expression: LD.FCS.Quadcopter.Throttle.bias
                                        * Referenced by: '<S223>/Bias'
                                        */
  real_T Integrator_IC_mg;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S220>/Integrator'
                                        */
  real_T Filter_IC_n;                  /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S220>/Filter'
                                        */
  real_T Integrator_IC_b2;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S216>/Integrator'
                                        */
  real_T Filter_IC_fw;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S216>/Filter'
                                        */
  real_T dYaw_Gain_p;                  /* Expression: LD.FCS.Quadcopter.Yaw.dYaw
                                        * Referenced by: '<S226>/dYaw'
                                        */
  real_T Saturation7_UpperSat_d;       /* Expression: 10
                                        * Referenced by: '<S166>/Saturation7'
                                        */
  real_T Saturation7_LowerSat_l;       /* Expression: -10
                                        * Referenced by: '<S166>/Saturation7'
                                        */
  real_T Integrator_IC_e1;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S204>/Integrator'
                                        */
  real_T Filter_IC_cy;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S204>/Filter'
                                        */
  real_T Integrator_IC_ob;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S218>/Integrator'
                                        */
  real_T Filter_IC_hf;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S218>/Filter'
                                        */
  real_T Integrator_IC_po;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S214>/Integrator'
                                        */
  real_T Filter_IC_oz;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S214>/Filter'
                                        */
  real_T dPitch_Gain_g;                /* Expression: LD.FCS.Quadcopter.Pitch.dPitch
                                        * Referenced by: '<S224>/dPitch'
                                        */
  real_T Saturation_UpperSat_m;        /* Expression: 1
                                        * Referenced by: '<S221>/Saturation'
                                        */
  real_T Saturation_LowerSat_p;        /* Expression: 0
                                        * Referenced by: '<S221>/Saturation'
                                        */
  real_T motor2position_Value_o[3];    /* Expression: LD.Propulsion.motor2.Position
                                        * Referenced by: '<S224>/motor2 position '
                                        */
  real_T Constant1_Value_bj[3];        /* Expression: LD.Inertia.CG
                                        * Referenced by: '<S224>/Constant1'
                                        */
  real_T motor2position1_Value_n[3];   /* Expression: LD.Propulsion.motor3.Position
                                        * Referenced by: '<S224>/motor2 position 1'
                                        */
  real_T Constant2_Value_c[3];         /* Expression: LD.Inertia.CG
                                        * Referenced by: '<S224>/Constant2'
                                        */
  real_T Gain_Gain_fs;                 /* Expression: -1
                                        * Referenced by: '<S226>/Gain'
                                        */
  real_T Gain1_Gain_o;                 /* Expression: -1
                                        * Referenced by: '<S225>/Gain1'
                                        */
  real_T Saturation1_UpperSat_m;       /* Expression: 1
                                        * Referenced by: '<S221>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_e;       /* Expression: 0
                                        * Referenced by: '<S221>/Saturation1'
                                        */
  real_T Saturation2_UpperSat_i;       /* Expression: 1
                                        * Referenced by: '<S221>/Saturation2'
                                        */
  real_T Saturation2_LowerSat_p;       /* Expression: 0
                                        * Referenced by: '<S221>/Saturation2'
                                        */
  real_T Gain1_Gain_f;                 /* Expression: -1
                                        * Referenced by: '<S226>/Gain1'
                                        */
  real_T Saturation3_UpperSat_i;       /* Expression: 1
                                        * Referenced by: '<S221>/Saturation3'
                                        */
  real_T Saturation3_LowerSat_f;       /* Expression: 0
                                        * Referenced by: '<S221>/Saturation3'
                                        */
  real_T Integrator_IC_ey;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S217>/Integrator'
                                        */
  real_T Filter_IC_od;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S217>/Filter'
                                        */
  real_T Filter_IC_cm1;                /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S203>/Filter'
                                        */
  real_T Integrator_IC_pf;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S203>/Integrator'
                                        */
  real_T Gain8_Gain_o;                 /* Expression: -1
                                        * Referenced by: '<S168>/Gain8'
                                        */
  real_T Gain_Gain_iw;                 /* Expression: pi/180
                                        * Referenced by: '<S168>/Gain'
                                        */
  real_T Gain1_Gain_kz;                /* Expression: pi/180
                                        * Referenced by: '<S168>/Gain1'
                                        */
  real_T Gain2_Gain_i;                 /* Expression: pi/180
                                        * Referenced by: '<S168>/Gain2'
                                        */
  real_T Gain3_Gain_c;                 /* Expression: pi/180
                                        * Referenced by: '<S168>/Gain3'
                                        */
  real_T Constant15_Value_j;           /* Expression: 20
                                        * Referenced by: '<S168>/Constant15'
                                        */
  real_T Saturation_UpperSat_o;        /* Expression: 25
                                        * Referenced by: '<S168>/Saturation'
                                        */
  real_T Saturation_LowerSat_j;        /* Expression: -25
                                        * Referenced by: '<S168>/Saturation'
                                        */
  real_T Integrator_IC_du;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S233>/Integrator'
                                        */
  real_T Filter_IC_eb;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S233>/Filter'
                                        */
  real_T Integrator_IC_p2w;            /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S235>/Integrator'
                                        */
  real_T Filter_IC_lf;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S235>/Filter'
                                        */
  real_T Integrator_IC_js;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S234>/Integrator'
                                        */
  real_T Filter_IC_iz;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S234>/Filter'
                                        */
  real_T Gain4_Gain_fu;                /* Expression: pi/180
                                        * Referenced by: '<S168>/Gain4'
                                        */
  real_T Constant14_Value_c;           /* Expression: 0
                                        * Referenced by: '<S168>/Constant14'
                                        */
  real_T Integrator_IC_bm;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S237>/Integrator'
                                        */
  real_T Filter_IC_ib;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S237>/Filter'
                                        */
  real_T Integrator_IC_ds;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S236>/Integrator'
                                        */
  real_T Filter_IC_ck;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S236>/Filter'
                                        */
  real_T Constant11_Value_d;           /* Expression: 0
                                        * Referenced by: '<S168>/Constant11'
                                        */
  real_T Gain5_Gain_k;                 /* Expression: pi/180
                                        * Referenced by: '<S168>/Gain5'
                                        */
  real_T Constant16_Value_m;           /* Expression: 0
                                        * Referenced by: '<S168>/Constant16'
                                        */
  real_T Integrator_IC_h;              /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S239>/Integrator'
                                        */
  real_T Filter_IC_oo;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S239>/Filter'
                                        */
  real_T Integrator_IC_ep;             /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S238>/Integrator'
                                        */
  real_T Filter_IC_pk;                 /* Expression: InitialConditionForFilter
                                        * Referenced by: '<S238>/Filter'
                                        */
  real_T Constant17_Value_b;           /* Expression: 20*pi/180
                                        * Referenced by: '<S168>/Constant17'
                                        */
  real_T Constant12_Value_g;           /* Expression: 20*pi/180
                                        * Referenced by: '<S168>/Constant12'
                                        */
  real_T Gain6_Gain_d;                 /* Expression: pi/180
                                        * Referenced by: '<S168>/Gain6'
                                        */
  real_T Constant13_Value_c;           /* Expression: 20*pi/180
                                        * Referenced by: '<S168>/Constant13'
                                        */
  real_T Gain7_Gain_e;                 /* Expression: pi/180
                                        * Referenced by: '<S168>/Gain7'
                                        */
  real_T Constant6_Value_n[2];         /* Expression: [0,0]
                                        * Referenced by: '<S168>/Constant6'
                                        */
  real_T Constant4_Value_g;            /* Expression: 1
                                        * Referenced by: '<S168>/Constant4'
                                        */
  real_T Constant5_Value_p;            /* Expression: 0
                                        * Referenced by: '<S168>/Constant5'
                                        */
  real_T Constant7_Value_c;            /* Expression: 0
                                        * Referenced by: '<S168>/Constant7'
                                        */
  real_T Constant8_Value_o;            /* Expression: 0
                                        * Referenced by: '<S168>/Constant8'
                                        */
  real_T Constant3_Value_l;            /* Expression: 0
                                        * Referenced by: '<S168>/Constant3'
                                        */
  real_T Bias_Bias_n;                  /* Expression: -90
                                        * Referenced by: '<S315>/Bias'
                                        */
  real_T Gain_Gain_l;                  /* Expression: -1
                                        * Referenced by: '<S315>/Gain'
                                        */
  real_T Bias1_Bias;                   /* Expression: +90
                                        * Referenced by: '<S315>/Bias1'
                                        */
  real_T Bias_Bias_g;                  /* Expression: 180
                                        * Referenced by: '<S318>/Bias'
                                        */
  real_T Bias1_Bias_a;                 /* Expression: -180
                                        * Referenced by: '<S318>/Bias1'
                                        */
  real_T Bias_Bias_e;                  /* Expression: 180
                                        * Referenced by: '<S316>/Bias'
                                        */
  real_T Bias1_Bias_ae;                /* Expression: -180
                                        * Referenced by: '<S316>/Bias1'
                                        */
  real_T Constant1_Value_k;            /* Expression: 0
                                        * Referenced by: '<S312>/Constant1'
                                        */
  real_T Constant_Value_pw;            /* Expression: 180
                                        * Referenced by: '<S312>/Constant'
                                        */
  real_T Bias_Bias_p;                  /* Expression: -90
                                        * Referenced by: '<S321>/Bias'
                                        */
  real_T Gain_Gain_dr;                 /* Expression: -1
                                        * Referenced by: '<S321>/Gain'
                                        */
  real_T Bias1_Bias_o;                 /* Expression: +90
                                        * Referenced by: '<S321>/Bias1'
                                        */
  real_T Constant2_Value_b;            /* Expression: 360
                                        * Referenced by: '<S324>/Constant2'
                                        */
  real_T Bias_Bias_g4;                 /* Expression: 180
                                        * Referenced by: '<S324>/Bias'
                                        */
  real_T Bias1_Bias_b;                 /* Expression: -180
                                        * Referenced by: '<S324>/Bias1'
                                        */
  real_T Constant2_Value_p;            /* Expression: 360
                                        * Referenced by: '<S322>/Constant2'
                                        */
  real_T Bias_Bias_i;                  /* Expression: 180
                                        * Referenced by: '<S322>/Bias'
                                        */
  real_T Bias1_Bias_d;                 /* Expression: -180
                                        * Referenced by: '<S322>/Bias1'
                                        */
  real_T VflightPrelookup_BreakpointsData[1000];/* Expression: LD.Propulsion.propeller1.V(1,:)
                                                 * Referenced by: '<S354>/Vflight Prelookup'
                                                 */
  real_T RPMPrelookup_BreakpointsData[15];/* Expression: LD.Propulsion.propeller1.RPM(:,1)
                                           * Referenced by: '<S354>/RPM Prelookup'
                                           */
  real_T TorqueInterpolationUsingPrelookup_Table[15000];/* Expression: LD.Propulsion.propeller1.Torque
                                                         * Referenced by: '<S354>/Torque Interpolation Using Prelookup'
                                                         */
  real_T VflightPrelookup_BreakpointsData_j[1000];/* Expression: LD.Propulsion.propeller2.V(1,:)
                                                   * Referenced by: '<S368>/Vflight Prelookup'
                                                   */
  real_T RPMPrelookup_BreakpointsData_n[12];/* Expression: LD.Propulsion.propeller2.RPM(:,1)
                                             * Referenced by: '<S368>/RPM Prelookup'
                                             */
  real_T TorqueInterpolationUsingPrelookup_Table_b[12000];/* Expression: LD.Propulsion.propeller2.Torque
                                                           * Referenced by: '<S368>/Torque Interpolation Using Prelookup'
                                                           */
  real_T VflightPrelookup_BreakpointsData_i[1000];/* Expression: LD.Propulsion.propeller3.V(1,:)
                                                   * Referenced by: '<S382>/Vflight Prelookup'
                                                   */
  real_T RPMPrelookup_BreakpointsData_f[12];/* Expression: LD.Propulsion.propeller3.RPM(:,1)
                                             * Referenced by: '<S382>/RPM Prelookup'
                                             */
  real_T TorqueInterpolationUsingPrelookup_Table_i[12000];/* Expression: LD.Propulsion.propeller3.Torque
                                                           * Referenced by: '<S382>/Torque Interpolation Using Prelookup'
                                                           */
  real_T VflightPrelookup_BreakpointsData_b[1000];/* Expression: LD.Propulsion.propeller4.V(1,:)
                                                   * Referenced by: '<S396>/Vflight Prelookup'
                                                   */
  real_T RPMPrelookup_BreakpointsData_m[12];/* Expression: LD.Propulsion.propeller4.RPM(:,1)
                                             * Referenced by: '<S396>/RPM Prelookup'
                                             */
  real_T TorqueInterpolationUsingPrelookup_Table_h[12000];/* Expression: LD.Propulsion.propeller4.Torque
                                                           * Referenced by: '<S396>/Torque Interpolation Using Prelookup'
                                                           */
  real_T VflightPrelookup_BreakpointsData_m[1000];/* Expression: LD.Propulsion.propeller5.V(1,:)
                                                   * Referenced by: '<S410>/Vflight Prelookup'
                                                   */
  real_T RPMPrelookup_BreakpointsData_ms[12];/* Expression: LD.Propulsion.propeller5.RPM(:,1)
                                              * Referenced by: '<S410>/RPM Prelookup'
                                              */
  real_T TorqueInterpolationUsingPrelookup_Table_p[12000];/* Expression: LD.Propulsion.propeller5.Torque
                                                           * Referenced by: '<S410>/Torque Interpolation Using Prelookup'
                                                           */
  real_T Bias_Bias_pa;                 /* Expression: -90
                                        * Referenced by: '<S507>/Bias'
                                        */
  real_T Gain_Gain_cn;                 /* Expression: -1
                                        * Referenced by: '<S507>/Gain'
                                        */
  real_T Bias1_Bias_m;                 /* Expression: +90
                                        * Referenced by: '<S507>/Bias1'
                                        */
  real_T Bias_Bias_k;                  /* Expression: 180
                                        * Referenced by: '<S510>/Bias'
                                        */
  real_T Bias1_Bias_do;                /* Expression: -180
                                        * Referenced by: '<S510>/Bias1'
                                        */
  real_T Bias_Bias_ok;                 /* Expression: 180
                                        * Referenced by: '<S508>/Bias'
                                        */
  real_T Bias1_Bias_i;                 /* Expression: -180
                                        * Referenced by: '<S508>/Bias1'
                                        */
  real_T Constant1_Value_e;            /* Expression: 0
                                        * Referenced by: '<S504>/Constant1'
                                        */
  real_T Constant_Value_h;             /* Expression: 180
                                        * Referenced by: '<S504>/Constant'
                                        */
  real_T Bias_Bias_pz;                 /* Expression: -90
                                        * Referenced by: '<S513>/Bias'
                                        */
  real_T Gain_Gain_bu;                 /* Expression: -1
                                        * Referenced by: '<S513>/Gain'
                                        */
  real_T Bias1_Bias_h;                 /* Expression: +90
                                        * Referenced by: '<S513>/Bias1'
                                        */
  real_T Constant2_Value_f;            /* Expression: 360
                                        * Referenced by: '<S516>/Constant2'
                                        */
  real_T Bias_Bias_i4;                 /* Expression: 180
                                        * Referenced by: '<S516>/Bias'
                                        */
  real_T Bias1_Bias_p;                 /* Expression: -180
                                        * Referenced by: '<S516>/Bias1'
                                        */
  real_T Constant2_Value_bb;           /* Expression: 360
                                        * Referenced by: '<S514>/Constant2'
                                        */
  real_T Bias_Bias_m;                  /* Expression: 180
                                        * Referenced by: '<S514>/Bias'
                                        */
  real_T Bias1_Bias_mf;                /* Expression: -180
                                        * Referenced by: '<S514>/Bias1'
                                        */
  real_T Trace3Phi0_Value[4];          /* Expression: [0 1 0 0]
                                        * Referenced by: '<S557>/Trace=3=>Phi=0'
                                        */
  real_T Shiftright_table[24];         /* Expression: [1 1 1;-1 1 1;1 1 -1;1 1 1;1 -1 1;1 1 1;1 1 1;1 1 1]
                                        * Referenced by: '<S559>/Shift right'
                                        */
  real_T Gain1_Gain_kf;                /* Expression: -1
                                        * Referenced by: '<S559>/Gain1'
                                        */
  real_T Constant_Value_hj;            /* Expression: 1
                                        * Referenced by: '<S558>/Constant'
                                        */
  real_T Gain_Gain_a;                  /* Expression: 0.5
                                        * Referenced by: '<S558>/Gain'
                                        */
  real_T Constant_Value_jg;            /* Expression: 0
                                        * Referenced by: '<S560>/Constant'
                                        */
  real_T Pi1_Value[3];                 /* Expression: [1 1 1]
                                        * Referenced by: '<S559>/Pi1'
                                        */
  real_T Switch_Threshold_p;           /* Expression: 0
                                        * Referenced by: '<S559>/Switch'
                                        */
  real_T Pi_Value;                     /* Expression: pi
                                        * Referenced by: '<S558>/Pi'
                                        */
  real_T Constant_Value_nm;            /* Expression: 1
                                        * Referenced by: '<S556>/Constant'
                                        */
  real_T Gain1_Gain_a;                 /* Expression: 0.5
                                        * Referenced by: '<S556>/Gain1'
                                        */
  real_T Gain_Gain_ao;                 /* Expression: 2
                                        * Referenced by: '<S556>/Gain'
                                        */
  real_T raddeg4_Gain;                 /* Expression: -1
                                        * Referenced by: '<S563>/rad-->deg4'
                                        */
  real_T IC4_Value_e[3];               /* Expression: [initialValues.Xe0 initialValues.Ye0 initialValues.Ze0]
                                        * Referenced by: '<S4>/IC4'
                                        */
  real_T Constant2_Value_l;            /* Expression: 1
                                        * Referenced by: '<S328>/Constant2'
                                        */
  real_T Constant1_Value_br;           /* Expression: R
                                        * Referenced by: '<S328>/Constant1'
                                        */
  real_T Constant_Value_c;             /* Expression: 1
                                        * Referenced by: '<S330>/Constant'
                                        */
  real_T Constant_Value_hu;            /* Expression: 1
                                        * Referenced by: '<S332>/Constant'
                                        */
  real_T Constant_Value_jd;            /* Expression: F
                                        * Referenced by: '<S331>/Constant'
                                        */
  real_T f_Value;                      /* Expression: 1
                                        * Referenced by: '<S331>/f'
                                        */
  real_T Constant_Value_p5;            /* Expression: 1
                                        * Referenced by: '<S328>/Constant'
                                        */
  real_T Constant3_Value_e;            /* Expression: 1
                                        * Referenced by: '<S328>/Constant3'
                                        */
  real_T Constant2_Value_g;            /* Expression: 360
                                        * Referenced by: '<S318>/Constant2'
                                        */
  real_T Constant_Value_f;             /* Expression: 180
                                        * Referenced by: '<S311>/Constant'
                                        */
  real_T Constant1_Value_a;            /* Expression: 0
                                        * Referenced by: '<S311>/Constant1'
                                        */
  real_T Constant2_Value_a;            /* Expression: 360
                                        * Referenced by: '<S316>/Constant2'
                                        */
  real_T Constant_Value_e;             /* Expression: initialValues.href
                                        * Referenced by: '<S4>/Constant'
                                        */
  real_T IC11_Value[3];                /* Expression: [initialValues.uned0 initialValues.vned0 initialValues.wned0]
                                        * Referenced by: '<S4>/IC11'
                                        */
  real_T IC12_Value[3];                /* Expression: [0 0 -initialValues.href]
                                        * Referenced by: '<S4>/IC12'
                                        */
  real_T IC2_Value_b[3];               /* Expression: [initialValues.roll0, initialValues.pitch0, initialValues.yaw0]
                                        * Referenced by: '<S4>/IC2'
                                        */
  real_T IC5_Value[3];                 /* Expression: [initialValues.u0 initialValues.v0 initialValues.w0]
                                        * Referenced by: '<S4>/IC5'
                                        */
  real_T IC3_Value_d[3];               /* Expression: [initialValues.p0 initialValues.q0 initialValues.r0]
                                        * Referenced by: '<S4>/IC3'
                                        */
  real_T TransferFcn_A;                /* Computed Parameter: TransferFcn_A
                                        * Referenced by: '<S250>/Transfer Fcn'
                                        */
  real_T TransferFcn_C;                /* Computed Parameter: TransferFcn_C
                                        * Referenced by: '<S250>/Transfer Fcn'
                                        */
  real_T TransferFcn1_A;               /* Computed Parameter: TransferFcn1_A
                                        * Referenced by: '<S250>/Transfer Fcn1'
                                        */
  real_T TransferFcn1_C;               /* Computed Parameter: TransferFcn1_C
                                        * Referenced by: '<S250>/Transfer Fcn1'
                                        */
  real_T TransferFcn2_A;               /* Computed Parameter: TransferFcn2_A
                                        * Referenced by: '<S250>/Transfer Fcn2'
                                        */
  real_T TransferFcn2_C;               /* Computed Parameter: TransferFcn2_C
                                        * Referenced by: '<S250>/Transfer Fcn2'
                                        */
  real_T TransferFcn3_A;               /* Computed Parameter: TransferFcn3_A
                                        * Referenced by: '<S250>/Transfer Fcn3'
                                        */
  real_T TransferFcn3_C;               /* Computed Parameter: TransferFcn3_C
                                        * Referenced by: '<S250>/Transfer Fcn3'
                                        */
  real_T TransferFcn4_A;               /* Computed Parameter: TransferFcn4_A
                                        * Referenced by: '<S250>/Transfer Fcn4'
                                        */
  real_T TransferFcn4_C;               /* Computed Parameter: TransferFcn4_C
                                        * Referenced by: '<S250>/Transfer Fcn4'
                                        */
  real_T TransferFcn5_A;               /* Computed Parameter: TransferFcn5_A
                                        * Referenced by: '<S250>/Transfer Fcn5'
                                        */
  real_T TransferFcn5_C;               /* Computed Parameter: TransferFcn5_C
                                        * Referenced by: '<S250>/Transfer Fcn5'
                                        */
  real_T TransferFcn1_A_l;             /* Computed Parameter: TransferFcn1_A_l
                                        * Referenced by: '<S249>/Transfer Fcn1'
                                        */
  real_T TransferFcn1_C_b;             /* Computed Parameter: TransferFcn1_C_b
                                        * Referenced by: '<S249>/Transfer Fcn1'
                                        */
  real_T TransferFcn4_A_f;             /* Computed Parameter: TransferFcn4_A_f
                                        * Referenced by: '<S249>/Transfer Fcn4'
                                        */
  real_T TransferFcn4_C_c;             /* Computed Parameter: TransferFcn4_C_c
                                        * Referenced by: '<S249>/Transfer Fcn4'
                                        */
  real_T TransferFcn5_A_m;             /* Computed Parameter: TransferFcn5_A_m
                                        * Referenced by: '<S249>/Transfer Fcn5'
                                        */
  real_T TransferFcn5_C_h;             /* Computed Parameter: TransferFcn5_C_h
                                        * Referenced by: '<S249>/Transfer Fcn5'
                                        */
  real_T IC7_Value;                    /* Expression: initialValues.alpha_wb
                                        * Referenced by: '<S4>/IC7'
                                        */
  real_T Constant_Value_ep;            /* Expression: 0
                                        * Referenced by: '<S245>/Constant'
                                        */
  real_T Switch_Threshold_o;           /* Expression: 0
                                        * Referenced by: '<S245>/Switch'
                                        */
  real_T IC8_Value;                    /* Expression: 0
                                        * Referenced by: '<S4>/IC8'
                                        */
  real_T RateLimiter_RisingLim;        /* Expression: 100
                                        * Referenced by: '<S4>/Rate Limiter'
                                        */
  real_T RateLimiter_FallingLim;       /* Expression: -100
                                        * Referenced by: '<S4>/Rate Limiter'
                                        */
  real_T RateLimiter1_RisingLim;       /* Expression: 100
                                        * Referenced by: '<S4>/Rate Limiter1'
                                        */
  real_T RateLimiter1_FallingLim;      /* Expression: -100
                                        * Referenced by: '<S4>/Rate Limiter1'
                                        */
  real_T masskg_Value;                 /* Expression: LD.Inertia.mass
                                        * Referenced by: '<S246>/mass (kg)'
                                        */
  real_T Integrator_IC_evh;            /* Expression: 0
                                        * Referenced by: '<S247>/Integrator'
                                        */
  real_T Constant_Value_co;            /* Expression: 3600
                                        * Referenced by: '<S247>/Constant'
                                        */
  real_T Battery_Capacity_Value;       /* Expression: 22
                                        * Referenced by: '<S343>/Battery_Capacity'
                                        */
  real_T profundidad_descarga_Value;   /* Expression: 0.95
                                        * Referenced by: '<S343>/profundidad_descarga'
                                        */
  real_T Gain_Gain_e;                  /* Expression: -1
                                        * Referenced by: '<S7>/Gain'
                                        */
  real_T Constant_Value_fu[2];         /* Expression: [maxzero 1]
                                        * Referenced by: '<S553>/Constant'
                                        */
  real_T Merge_InitialOutput_k;        /* Computed Parameter: Merge_InitialOutput_k
                                        * Referenced by: '<S553>/Merge'
                                        */
  real_T raddeg_Gain;                  /* Expression: -1
                                        * Referenced by: '<S551>/rad-->deg'
                                        */
  real_T Constant_Value_jm;            /* Expression: 1.01
                                        * Referenced by: '<S551>/Constant'
                                        */
  real_T raddeg1_Gain;                 /* Expression: -1
                                        * Referenced by: '<S551>/rad-->deg1'
                                        */
  real_T Constant_Value_dz;            /* Expression: eps
                                        * Referenced by: '<S563>/Constant'
                                        */
  real_T Switch_Threshold_m;           /* Expression: 0
                                        * Referenced by: '<S563>/Switch'
                                        */
  real_T InitialJuliandate_Value;      /* Expression: juliandate(datetime(initialValues.Date,'InputFormat','dd-MMM-yyyy hh:mm:ss'))
                                        * Referenced by: '<S2>/Initial Julian  date'
                                        */
  real_T kclock_Gain;                  /* Expression: 1/(60*60*24*(365 + leapyear(floor(decyear(initialValues.Date,'dd-mmm-yyyy')))))
                                        * Referenced by: '<S2>/kclock'
                                        */
  real_T GravityinEarthAxes_Gain[3];   /* Expression: [0 0 LD.Inertia.mass]
                                        * Referenced by: '<S2>/Gravity in Earth Axes'
                                        */
  real_T InitialDecimaldate_Value;     /* Expression: decyear(datetime(initialValues.Date,'InputFormat','dd-MMM-yyyy hh:mm:ss'))
                                        * Referenced by: '<S2>/Initial Decimal date'
                                        */
  real_T motor2position_Value_c[3];    /* Expression: LD.LandingGear.A.Position
                                        * Referenced by: '<S13>/motor2 position '
                                        */
  real_T Constant4_Value_p;            /* Expression: 0
                                        * Referenced by: '<S28>/Constant4'
                                        */
  real_T Constant5_Value_l;            /* Expression: 0
                                        * Referenced by: '<S28>/Constant5'
                                        */
  real_T Constant3_Value_lv;           /* Expression: 0
                                        * Referenced by: '<S28>/Constant3'
                                        */
  real_T hground_value1_Value;         /* Expression: initialValues.href
                                        * Referenced by: '<S28>/hground_value1'
                                        */
  real_T Switch3_Threshold;            /* Expression: 0
                                        * Referenced by: '<S28>/Switch3'
                                        */
  real_T motor2position1_Value_m[3];   /* Expression: LD.LandingGear.B.Position
                                        * Referenced by: '<S13>/motor2 position 1'
                                        */
  real_T Constant4_Value_c;            /* Expression: 0
                                        * Referenced by: '<S29>/Constant4'
                                        */
  real_T Constant5_Value_b;            /* Expression: 0
                                        * Referenced by: '<S29>/Constant5'
                                        */
  real_T Constant3_Value_k;            /* Expression: 0
                                        * Referenced by: '<S29>/Constant3'
                                        */
  real_T hground_value1_Value_o;       /* Expression: initialValues.href
                                        * Referenced by: '<S29>/hground_value1'
                                        */
  real_T Switch3_Threshold_g;          /* Expression: 0
                                        * Referenced by: '<S29>/Switch3'
                                        */
  real_T motor2position2_Value[3];     /* Expression: LD.LandingGear.C.Position
                                        * Referenced by: '<S13>/motor2 position 2'
                                        */
  real_T Constant4_Value_c3;           /* Expression: 0
                                        * Referenced by: '<S30>/Constant4'
                                        */
  real_T Constant5_Value_d;            /* Expression: 0
                                        * Referenced by: '<S30>/Constant5'
                                        */
  real_T Constant3_Value_h;            /* Expression: 0
                                        * Referenced by: '<S30>/Constant3'
                                        */
  real_T hground_value1_Value_h;       /* Expression: initialValues.href
                                        * Referenced by: '<S30>/hground_value1'
                                        */
  real_T Switch3_Threshold_b;          /* Expression: 0
                                        * Referenced by: '<S30>/Switch3'
                                        */
  real_T IC_Value_e[3];                /* Expression: [0 0 0]
                                        * Referenced by: '<S13>/IC'
                                        */
  real_T IC1_Value_b[3];               /* Expression: [0 0 0]
                                        * Referenced by: '<S13>/IC1'
                                        */
  real_T IC2_Value_c[3];               /* Expression: [0 0 0]
                                        * Referenced by: '<S13>/IC2'
                                        */
  real_T IC3_Value_k[3];               /* Expression: [0 0 0]
                                        * Referenced by: '<S13>/IC3'
                                        */
  real_T IC4_Value_g[3];               /* Expression: [0 0 0]
                                        * Referenced by: '<S13>/IC4'
                                        */
  real_T IC5_Value_b[3];               /* Expression: [0 0 0]
                                        * Referenced by: '<S13>/IC5'
                                        */
  real_T epoch_Value;                  /* Expression: epoch
                                        * Referenced by: '<S104>/epoch'
                                        */
  real_T otime_X0;                     /* Expression: -1000
                                        * Referenced by: '<S112>/otime'
                                        */
  real_T u80deg_UpperSat;              /* Expression: 180.0
                                        * Referenced by: '<S16>/+//- 180 deg'
                                        */
  real_T u80deg_LowerSat;              /* Expression: -180.0
                                        * Referenced by: '<S16>/+//- 180 deg'
                                        */
  real_T u0deg_UpperSat;               /* Expression: 90.0
                                        * Referenced by: '<S16>/+//- 90 deg'
                                        */
  real_T u0deg_LowerSat;               /* Expression: -90.0
                                        * Referenced by: '<S16>/+//- 90 deg'
                                        */
  real_T olon_X0;                      /* Expression: -1000
                                        * Referenced by: '<S111>/olon'
                                        */
  real_T olat_X0;                      /* Expression: -1000
                                        * Referenced by: '<S110>/olat'
                                        */
  real_T uto1000000m_UpperSat;         /* Expression: 1000000.0
                                        * Referenced by: '<S16>/0 to 1,000,000 m'
                                        */
  real_T uto1000000m_LowerSat;         /* Expression: 0
                                        * Referenced by: '<S16>/0 to 1,000,000 m'
                                        */
  real_T Gain_Gain_il;                 /* Expression: 0.001
                                        * Referenced by: '<S16>/Gain'
                                        */
  real_T oalt_X0;                      /* Expression: -1000
                                        * Referenced by: '<S110>/oalt'
                                        */
  real_T re_Value;                     /* Expression: 6371.2
                                        * Referenced by: '<S104>/re'
                                        */
  real_T Constant_Value_bp[3];         /* Expression: [0,0,0]
                                        * Referenced by: '<S15>/Constant'
                                        */
  real_T Constant2_Value_cd[3];        /* Expression: [0,0,0]
                                        * Referenced by: '<S15>/Constant2'
                                        */
  real_T u_Value;                      /* Expression: 1.0
                                        * Referenced by: '<S47>/2'
                                        */
  real_T LimitFunction10ftto1000ft_UpperSat;/* Expression: max_height_low
                                             * Referenced by: '<S91>/Limit Function 10ft to 1000ft'
                                             */
  real_T LimitFunction10ftto1000ft_LowerSat;/* Expression: 10
                                             * Referenced by: '<S91>/Limit Function 10ft to 1000ft'
                                             */
  real_T Lw_Gain;                      /* Expression: 1
                                        * Referenced by: '<S62>/Lw'
                                        */
  real_T sigma_wg_Gain;                /* Expression: 0.1
                                        * Referenced by: '<S74>/sigma_wg '
                                        */
  real_T PreLookUpIndexSearchaltitude_BreakpointsData[12];/* Expression: h_vec
                                                           * Referenced by: '<S73>/PreLook-Up Index Search  (altitude)'
                                                           */
  real_T PreLookUpIndexSearchprobofexceed_BreakpointsData[7];/* Expression: [1:7]
                                                              * Referenced by: '<S73>/PreLook-Up Index Search  (prob of exceed)'
                                                              */
  real_T MediumHighAltitudeIntensity_Table[84];/* Expression: sigma_vec'
                                                * Referenced by: '<S73>/Medium//High Altitude Intensity'
                                                */
  real_T WhiteNoise_Mean;              /* Expression: 0
                                        * Referenced by: '<S65>/White Noise'
                                        */
  real_T WhiteNoise_StdDev;            /* Computed Parameter: WhiteNoise_StdDev
                                        * Referenced by: '<S65>/White Noise'
                                        */
  real_T WhiteNoise_Seed;              /* Expression: seed
                                        * Referenced by: '<S65>/White Noise'
                                        */
  real_T Output_Gain;                  /* Expression: [sqrt(Cov)]/[sqrt(Ts)]
                                        * Referenced by: '<S65>/Output'
                                        */
  real_T WhiteNoise_Mean_o;            /* Expression: 0
                                        * Referenced by: '<S66>/White Noise'
                                        */
  real_T WhiteNoise_StdDev_c;          /* Computed Parameter: WhiteNoise_StdDev_c
                                        * Referenced by: '<S66>/White Noise'
                                        */
  real_T WhiteNoise_Seed_f[3];         /* Expression: seed
                                        * Referenced by: '<S66>/White Noise'
                                        */
  real_T Output_Gain_g[3];             /* Expression: [sqrt(Cov)]/[sqrt(Ts)]
                                        * Referenced by: '<S66>/Output'
                                        */
  real_T LimitHeighth1000ft_UpperSat;  /* Expression: max_height_low
                                        * Referenced by: '<S74>/Limit Height h<1000ft'
                                        */
  real_T LimitHeighth1000ft_LowerSat;  /* Expression: 0
                                        * Referenced by: '<S74>/Limit Height h<1000ft'
                                        */
  real_T Lv_Gain;                      /* Expression: 1
                                        * Referenced by: '<S62>/Lv'
                                        */
  real_T ref_heightz0_Value;           /* Expression: 20/z0
                                        * Referenced by: '<S49>/ref_height//z0'
                                        */
  real_T Wdeg1_Value;                  /* Expression: 0
                                        * Referenced by: '<S49>/Wdeg1'
                                        */
  real_T Constant1_Value_e2;           /* Expression: 0
                                        * Referenced by: '<S3>/Constant1'
                                        */
  real_T ComandedpositioninEarthaxis_Value[3];/* Expression: [10,2,20]
                                               * Referenced by: '<S5>/Comanded position in Earth axis'
                                               */
  real_T JoystickInput_P1_Size[2];     /* Computed Parameter: JoystickInput_P1_Size
                                        * Referenced by: '<S5>/Joystick Input'
                                        */
  real_T JoystickInput_P1;             /* Expression: joyid
                                        * Referenced by: '<S5>/Joystick Input'
                                        */
  real_T JoystickInput_P2_Size[2];     /* Computed Parameter: JoystickInput_P2_Size
                                        * Referenced by: '<S5>/Joystick Input'
                                        */
  real_T JoystickInput_P2;             /* Expression: adjustports
                                        * Referenced by: '<S5>/Joystick Input'
                                        */
  real_T JoystickInput_P3_Size[2];     /* Computed Parameter: JoystickInput_P3_Size
                                        * Referenced by: '<S5>/Joystick Input'
                                        */
  real_T JoystickInput_P3;             /* Expression: forcefeed
                                        * Referenced by: '<S5>/Joystick Input'
                                        */
  real_T Constant2_Value_lr;           /* Expression: 360
                                        * Referenced by: '<S510>/Constant2'
                                        */
  real_T Constant_Value_pu;            /* Expression: 180
                                        * Referenced by: '<S503>/Constant'
                                        */
  real_T Constant1_Value_g;            /* Expression: 0
                                        * Referenced by: '<S503>/Constant1'
                                        */
  real_T Constant2_Value_cg;           /* Expression: 360
                                        * Referenced by: '<S508>/Constant2'
                                        */
  real_T Constant2_Value_n;            /* Expression: 1
                                        * Referenced by: '<S520>/Constant2'
                                        */
  real_T Constant1_Value_m;            /* Expression: R
                                        * Referenced by: '<S520>/Constant1'
                                        */
  real_T Constant_Value_ke;            /* Expression: 1
                                        * Referenced by: '<S522>/Constant'
                                        */
  real_T Constant_Value_bf;            /* Expression: 1
                                        * Referenced by: '<S524>/Constant'
                                        */
  real_T Constant_Value_gt;            /* Expression: F
                                        * Referenced by: '<S523>/Constant'
                                        */
  real_T f_Value_n;                    /* Expression: 1
                                        * Referenced by: '<S523>/f'
                                        */
  real_T Constant_Value_pj;            /* Expression: 1
                                        * Referenced by: '<S520>/Constant'
                                        */
  real_T Constant3_Value_f;            /* Expression: 1
                                        * Referenced by: '<S520>/Constant3'
                                        */
  real_T Constant13_Value_i;           /* Expression: initialValues.href
                                        * Referenced by: '<S498>/Constant13'
                                        */
  real_T TransferFcnX_A[2];            /* Computed Parameter: TransferFcnX_A
                                        * Referenced by: '<S538>/Transfer Fcn X'
                                        */
  real_T TransferFcnX_C[2];            /* Computed Parameter: TransferFcnX_C
                                        * Referenced by: '<S538>/Transfer Fcn X'
                                        */
  real_T TransferFcnY_A[2];            /* Computed Parameter: TransferFcnY_A
                                        * Referenced by: '<S538>/Transfer Fcn Y'
                                        */
  real_T TransferFcnY_C[2];            /* Computed Parameter: TransferFcnY_C
                                        * Referenced by: '<S538>/Transfer Fcn Y'
                                        */
  real_T TransferFcnZ_A[2];            /* Computed Parameter: TransferFcnZ_A
                                        * Referenced by: '<S538>/Transfer Fcn Z'
                                        */
  real_T TransferFcnZ_C[2];            /* Computed Parameter: TransferFcnZ_C
                                        * Referenced by: '<S538>/Transfer Fcn Z'
                                        */
  real_T Constant_Value_ob;            /* Expression: dtype_a
                                        * Referenced by: '<S533>/Constant'
                                        */
  real_T ZeroOrderHold1_Gain;          /* Expression: 1
                                        * Referenced by: '<S531>/Zero-Order Hold1'
                                        */
  real_T ZeroOrderHold2_Gain;          /* Expression: 1
                                        * Referenced by: '<S531>/Zero-Order Hold2'
                                        */
  real_T ZeroOrderHold_Gain;           /* Expression: 1
                                        * Referenced by: '<S531>/Zero-Order Hold'
                                        */
  real_T ZeroOrderHold4_Gain;          /* Expression: 1
                                        * Referenced by: '<S531>/Zero-Order Hold4'
                                        */
  real_T Gain_Gain_iwz[3];             /* Expression: [1 -1 1]
                                        * Referenced by: '<S531>/Gain'
                                        */
  real_T ZeroOrderHold3_Gain;          /* Expression: 1
                                        * Referenced by: '<S531>/Zero-Order Hold3'
                                        */
  real_T Switch_Threshold_o5;          /* Expression: 0.5
                                        * Referenced by: '<S533>/Switch'
                                        */
  real_T WhiteNoise_Mean_d;            /* Expression: 0
                                        * Referenced by: '<S534>/White Noise'
                                        */
  real_T WhiteNoise_StdDev_n;          /* Computed Parameter: WhiteNoise_StdDev_n
                                        * Referenced by: '<S534>/White Noise'
                                        */
  real_T WhiteNoise_Seed_n[3];         /* Expression: seed
                                        * Referenced by: '<S534>/White Noise'
                                        */
  real_T Output_Gain_o[3];             /* Expression: [sqrt(Cov)]/[sqrt(Ts)]
                                        * Referenced by: '<S534>/Output'
                                        */
  real_T Saturation_UpperSat_a[3];     /* Expression: a_sath
                                        * Referenced by: '<S531>/Saturation'
                                        */
  real_T Saturation_LowerSat_o[3];     /* Expression: a_satl
                                        * Referenced by: '<S531>/Saturation'
                                        */
  real_T TransferFcnX_A_d[2];          /* Computed Parameter: TransferFcnX_A_d
                                        * Referenced by: '<S550>/Transfer Fcn X'
                                        */
  real_T TransferFcnX_C_d[2];          /* Computed Parameter: TransferFcnX_C_d
                                        * Referenced by: '<S550>/Transfer Fcn X'
                                        */
  real_T TransferFcnY_A_h[2];          /* Computed Parameter: TransferFcnY_A_h
                                        * Referenced by: '<S550>/Transfer Fcn Y'
                                        */
  real_T TransferFcnY_C_p[2];          /* Computed Parameter: TransferFcnY_C_p
                                        * Referenced by: '<S550>/Transfer Fcn Y'
                                        */
  real_T TransferFcnZ_A_p[2];          /* Computed Parameter: TransferFcnZ_A_p
                                        * Referenced by: '<S550>/Transfer Fcn Z'
                                        */
  real_T TransferFcnZ_C_d[2];          /* Computed Parameter: TransferFcnZ_C_d
                                        * Referenced by: '<S550>/Transfer Fcn Z'
                                        */
  real_T Constant_Value_id;            /* Expression: dtype_g
                                        * Referenced by: '<S547>/Constant'
                                        */
  real_T ZeroOrderHold_Gain_c;         /* Expression: 1
                                        * Referenced by: '<S532>/Zero-Order Hold'
                                        */
  real_T ZeroOrderHold1_Gain_k;        /* Expression: 1
                                        * Referenced by: '<S532>/Zero-Order Hold1'
                                        */
  real_T Switch_Threshold_f;           /* Expression: 0.5
                                        * Referenced by: '<S547>/Switch'
                                        */
  real_T WhiteNoise_Mean_n;            /* Expression: 0
                                        * Referenced by: '<S548>/White Noise'
                                        */
  real_T WhiteNoise_StdDev_l;          /* Computed Parameter: WhiteNoise_StdDev_l
                                        * Referenced by: '<S548>/White Noise'
                                        */
  real_T WhiteNoise_Seed_a[3];         /* Expression: seed
                                        * Referenced by: '<S548>/White Noise'
                                        */
  real_T Output_Gain_h[3];             /* Expression: [sqrt(Cov)]/[sqrt(Ts)]
                                        * Referenced by: '<S548>/Output'
                                        */
  real_T Saturation_UpperSat_c[3];     /* Expression: g_sath
                                        * Referenced by: '<S532>/Saturation'
                                        */
  real_T Saturation_LowerSat_m[3];     /* Expression: g_satl
                                        * Referenced by: '<S532>/Saturation'
                                        */
  real_T phithetapsi_IC[3];            /* Expression: [initialValues.roll0 initialValues.pitch0 initialValues.yaw0]
                                        * Referenced by: '<S525>/phi theta psi'
                                        */
  real_T Constant1_Value_i;            /* Expression: 0
                                        * Referenced by: '<S259>/Constant1'
                                        */
  real_T Constant_Value_bw;            /* Expression: 0
                                        * Referenced by: '<S259>/Constant'
                                        */
  real_T Constant2_Value_ne;           /* Expression: 1
                                        * Referenced by: '<S259>/Constant2'
                                        */
  real_T Constant4_Value_p1;           /* Expression: 0
                                        * Referenced by: '<S259>/Constant4'
                                        */
  real_T Constant3_Value_p;            /* Expression: 0
                                        * Referenced by: '<S259>/Constant3'
                                        */
  real_T preAlpha1_BreakpointsData[9]; /* Expression: LD.alpha
                                        * Referenced by: '<S242>/preAlpha1'
                                        */
  real_T LDflagsalphalength_Value;     /* Expression: LD.flags.alpha.length
                                        * Referenced by: '<S242>/LD.flags.alpha.length'
                                        */
  real_T Constant1_Value_dt;           /* Expression: 0
                                        * Referenced by: '<S242>/Constant 1'
                                        */
  real_T Switch7_Threshold;            /* Expression: checks.inconsistentLength
                                        * Referenced by: '<S242>/Switch7'
                                        */
  real_T preBeta_BreakpointsData[2];   /* Expression: LD.beta
                                        * Referenced by: '<S242>/preBeta'
                                        */
  real_T LDflagsbetalength_Value;      /* Expression: LD.flags.beta.length
                                        * Referenced by: '<S242>/LD.flags.beta.length'
                                        */
  real_T Constant15_Value_js;          /* Expression: 0
                                        * Referenced by: '<S242>/Constant 15'
                                        */
  real_T Switch1_Threshold_f;          /* Expression: checks.inconsistentLength
                                        * Referenced by: '<S242>/Switch1'
                                        */
  real_T hlookup1_BreakpointsData[2];  /* Expression: LD.alt
                                        * Referenced by: '<S242>/h  look-up1'
                                        */
  real_T LDflagsaltlength_Value;       /* Expression: LD.flags.alt.length
                                        * Referenced by: '<S242>/LD.flags.alt.length'
                                        */
  real_T Constant3_Value_b;            /* Expression: 0
                                        * Referenced by: '<S242>/Constant 3'
                                        */
  real_T Switch8_Threshold;            /* Expression: checks.inconsistentLength
                                        * Referenced by: '<S242>/Switch8'
                                        */
  real_T lookup1_BreakpointsData[2];   /* Expression: LD.xcg
                                        * Referenced by: '<S242>/ look-up1'
                                        */
  real_T LDflagsxcglength_Value;       /* Expression: LD.flags.xcg.length
                                        * Referenced by: '<S242>/LD.flags.xcg.length'
                                        */
  real_T Constant5_Value_n;            /* Expression: 0
                                        * Referenced by: '<S242>/Constant 5'
                                        */
  real_T Switch9_Threshold;            /* Expression: checks.inconsistentLength
                                        * Referenced by: '<S242>/Switch9'
                                        */
  real_T p1_BreakpointsData[2];        /* Expression: LD.deltae
                                        * Referenced by: '<S242>/p1'
                                        */
  real_T LDflagsdeltaelength_Value;    /* Expression: LD.flags.deltae.length
                                        * Referenced by: '<S242>/LD.flags.deltae.length'
                                        */
  real_T Constant7_Value_o;            /* Expression: 0
                                        * Referenced by: '<S242>/Constant 7'
                                        */
  real_T Switch10_Threshold;           /* Expression: checks.inconsistentLength
                                        * Referenced by: '<S242>/Switch10'
                                        */
  real_T lookup1_BreakpointsData_g[2]; /* Expression: LD.deltar
                                        * Referenced by: '<S242>/look-up1'
                                        */
  real_T LDflagsdeltarlength_Value;    /* Expression: LD.flags.deltar.length
                                        * Referenced by: '<S242>/LD.flags.deltar.length'
                                        */
  real_T Constant9_Value_a;            /* Expression: 0
                                        * Referenced by: '<S242>/Constant 9'
                                        */
  real_T Switch11_Threshold;           /* Expression: checks.inconsistentLength
                                        * Referenced by: '<S242>/Switch11'
                                        */
  real_T frlookup1_BreakpointsData[2]; /* Expression: LD.deltafr
                                        * Referenced by: '<S242>/fr look-up1'
                                        */
  real_T LDflagsdeltafrlength_Value;   /* Expression: LD.flags.deltafr.length
                                        * Referenced by: '<S242>/LD.flags.deltafr.length'
                                        */
  real_T Constant11_Value_dd;          /* Expression: 0
                                        * Referenced by: '<S242>/Constant 11'
                                        */
  real_T Switch12_Threshold;           /* Expression: checks.inconsistentLength
                                        * Referenced by: '<S242>/Switch12'
                                        */
  real_T prelookups1_BreakpointsData[2];/* Expression: LD.deltafl
                                         * Referenced by: '<S242>/pre look-ups1'
                                         */
  real_T LDflagsdeltafllength_Value;   /* Expression: LD.flags.deltafl.length
                                        * Referenced by: '<S242>/LD.flags.deltafl.length'
                                        */
  real_T Constant13_Value_cc;          /* Expression: 0
                                        * Referenced by: '<S242>/Constant 13'
                                        */
  real_T Switch13_Threshold;           /* Expression: checks.inconsistentLength
                                        * Referenced by: '<S242>/Switch13'
                                        */
  real_T CD0_Table[1152];              /* Expression: LD.Stability.CD0
                                        * Referenced by: '<S253>/CD0'
                                        */
  real_T CDalpha_Table[1152];          /* Expression: LD.Stability.CDalpha
                                        * Referenced by: '<S253>/CDalpha'
                                        */
  real_T CDalpha_dot_Table[1152];      /* Expression: LD.Stability.CDalpha_dot
                                        * Referenced by: '<S253>/CDalpha_dot'
                                        */
  real_T CDq_Table[1152];              /* Expression: LD.Stability.CDq
                                        * Referenced by: '<S253>/CDq'
                                        */
  real_T CDdeltae_Table[1152];         /* Expression: LD.Stability.CDdeltae
                                        * Referenced by: '<S253>/CDdeltae'
                                        */
  real_T CDdeltafr_Table[1152];        /* Expression: LD.Stability.CDdeltafr
                                        * Referenced by: '<S253>/CDdeltafr'
                                        */
  real_T CDdeltafrl_Table[1152];       /* Expression: LD.Stability.CDdeltafl
                                        * Referenced by: '<S253>/CDdeltafrl'
                                        */
  real_T CYbeta_Table[1152];           /* Expression: LD.Stability.CYbeta
                                        * Referenced by: '<S255>/CYbeta'
                                        */
  real_T CYp_Table[1152];              /* Expression: LD.Stability.CYp
                                        * Referenced by: '<S255>/CYp'
                                        */
  real_T CYr_Table[1152];              /* Expression: LD.Stability.CYr
                                        * Referenced by: '<S255>/CYr'
                                        */
  real_T CYdeltar_Table[1152];         /* Expression: LD.Stability.CYdeltar
                                        * Referenced by: '<S255>/CYdeltar '
                                        */
  real_T CYdeltafr_Table[1152];        /* Expression: LD.Stability.CYdeltafr
                                        * Referenced by: '<S255>/CYdeltafr'
                                        */
  real_T CYdeltafl_Table[1152];        /* Expression: LD.Stability.CYdeltafl
                                        * Referenced by: '<S255>/CYdeltafl'
                                        */
  real_T CL0_Table[1152];              /* Expression: LD.Stability.CL0
                                        * Referenced by: '<S254>/CL0'
                                        */
  real_T CLalpha_Table[1152];          /* Expression: LD.Stability.CLalpha
                                        * Referenced by: '<S254>/CLalpha'
                                        */
  real_T CLalpha_dot_Table[1152];      /* Expression: LD.Stability.CLalpha_dot
                                        * Referenced by: '<S254>/CLalpha_dot'
                                        */
  real_T CLq_Table[1152];              /* Expression: LD.Stability.CLq
                                        * Referenced by: '<S254>/CLq'
                                        */
  real_T CLdeltae_Table[1152];         /* Expression: LD.Stability.CLdeltae
                                        * Referenced by: '<S254>/CLdeltae'
                                        */
  real_T CLdeltafr_Table[1152];        /* Expression: LD.Stability.CLdeltafr
                                        * Referenced by: '<S254>/CLdeltafr'
                                        */
  real_T CLdeltafrl_Table[1152];       /* Expression: LD.Stability.CLdeltafl
                                        * Referenced by: '<S254>/CLdeltafrl'
                                        */
  real_T u2rhoV2_Gain;                 /* Expression: 1/2
                                        * Referenced by: '<S260>/1//2rhoV^2'
                                        */
  real_T coefAdjust_Gain[3];           /* Expression: [-1  1 -1]
                                        * Referenced by: '<S252>/coefAdjust'
                                        */
  real_T Constant1_Value_f3;           /* Expression: 0
                                        * Referenced by: '<S268>/Constant1'
                                        */
  real_T Constant_Value_bv;            /* Expression: 0
                                        * Referenced by: '<S268>/Constant'
                                        */
  real_T Constant2_Value_lrz;          /* Expression: 1
                                        * Referenced by: '<S268>/Constant2'
                                        */
  real_T Constant4_Value_i;            /* Expression: 0
                                        * Referenced by: '<S268>/Constant4'
                                        */
  real_T Constant3_Value_o;            /* Expression: 0
                                        * Referenced by: '<S268>/Constant3'
                                        */
  real_T Constant_Value_e5;            /* Expression: 0
                                        * Referenced by: '<S274>/Constant'
                                        */
  real_T Constant1_Value_lq;           /* Expression: 0
                                        * Referenced by: '<S274>/Constant1'
                                        */
  real_T Constant2_Value_nf;           /* Expression: 1
                                        * Referenced by: '<S274>/Constant2'
                                        */
  real_T Constant3_Value_mn;           /* Expression: 0
                                        * Referenced by: '<S274>/Constant3'
                                        */
  real_T Constant4_Value_pb;           /* Expression: 0
                                        * Referenced by: '<S274>/Constant4'
                                        */
  real_T Constant_Value_dh;            /* Expression: 0
                                        * Referenced by: '<S280>/Constant'
                                        */
  real_T Constant1_Value_j;            /* Expression: 0
                                        * Referenced by: '<S280>/Constant1'
                                        */
  real_T Constant2_Value_ho;           /* Expression: 1
                                        * Referenced by: '<S280>/Constant2'
                                        */
  real_T Constant3_Value_i;            /* Expression: 0
                                        * Referenced by: '<S280>/Constant3'
                                        */
  real_T Constant4_Value_k;            /* Expression: 0
                                        * Referenced by: '<S280>/Constant4'
                                        */
  real_T Clbeta_Table[1152];           /* Expression: LD.Stability.Clbeta
                                        * Referenced by: '<S256>/Clbeta'
                                        */
  real_T Clp_Table[1152];              /* Expression: LD.Stability.Clp
                                        * Referenced by: '<S256>/Clp'
                                        */
  real_T Clr_Table[1152];              /* Expression: LD.Stability.Clr
                                        * Referenced by: '<S256>/Clr'
                                        */
  real_T Cldeltar_Table[1152];         /* Expression: LD.Stability.Cldeltar
                                        * Referenced by: '<S256>/Cldeltar '
                                        */
  real_T Cldeltafr_Table[1152];        /* Expression: LD.Stability.Cldeltafr
                                        * Referenced by: '<S256>/Cldeltafr'
                                        */
  real_T Cldeltafl_Table[1152];        /* Expression: LD.Stability.Cldeltafl
                                        * Referenced by: '<S256>/Cldeltafl'
                                        */
  real_T Cm0_Table[1152];              /* Expression: LD.Stability.Cm0
                                        * Referenced by: '<S257>/Cm0'
                                        */
  real_T Cmalpha_Table[1152];          /* Expression: LD.Stability.Cmalpha
                                        * Referenced by: '<S257>/Cmalpha'
                                        */
  real_T Cmalpha_dot_Table[1152];      /* Expression: LD.Stability.Cmalpha_dot
                                        * Referenced by: '<S257>/Cmalpha_dot'
                                        */
  real_T Cmq_Table[1152];              /* Expression: LD.Stability.Cmq
                                        * Referenced by: '<S257>/Cmq'
                                        */
  real_T Cmdeltae_Table[1152];         /* Expression: LD.Stability.Cmdeltae
                                        * Referenced by: '<S257>/Cmdeltae'
                                        */
  real_T Cmdeltafr_Table[1152];        /* Expression: LD.Stability.Cmdeltafr
                                        * Referenced by: '<S257>/Cmdeltafr'
                                        */
  real_T Cmdeltafl_Table[1152];        /* Expression: LD.Stability.CLdeltafl
                                        * Referenced by: '<S257>/Cmdeltafl'
                                        */
  real_T Cnbeta_Table[1152];           /* Expression: LD.Stability.Cnbeta
                                        * Referenced by: '<S258>/Cnbeta'
                                        */
  real_T Cnp_Table[1152];              /* Expression: LD.Stability.Cnp
                                        * Referenced by: '<S258>/Cnp'
                                        */
  real_T Cnr_Table[1152];              /* Expression: LD.Stability.Cnr
                                        * Referenced by: '<S258>/Cnr'
                                        */
  real_T Cndeltar_Table[1152];         /* Expression: LD.Stability.Cndeltar
                                        * Referenced by: '<S258>/Cndeltar '
                                        */
  real_T Cndeltafr_Table[1152];        /* Expression: LD.Stability.Cndeltafr
                                        * Referenced by: '<S258>/Cndeltafr'
                                        */
  real_T Cndeltafl_Table[1152];        /* Expression: LD.Stability.Cndeltafl
                                        * Referenced by: '<S258>/Cndeltafl'
                                        */
  real_T IC1_Value_d[3];               /* Expression: [0 0 0]
                                        * Referenced by: '<S242>/IC1'
                                        */
  real_T IC2_Value_a[3];               /* Expression: [0 0 0]
                                        * Referenced by: '<S242>/IC2'
                                        */
  real_T Ixx_Value;                    /* Expression: LD.Inertia.Ixx
                                        * Referenced by: '<S246>/Ixx'
                                        */
  real_T Ixy_Value;                    /* Expression: LD.Inertia.Ixy
                                        * Referenced by: '<S246>/Ixy'
                                        */
  real_T Ixz_Value;                    /* Expression: LD.Inertia.Ixz
                                        * Referenced by: '<S246>/Ixz'
                                        */
  real_T Iyy_Value;                    /* Expression: LD.Inertia.Iyy
                                        * Referenced by: '<S246>/Iyy'
                                        */
  real_T Iyz_Value;                    /* Expression: LD.Inertia.Iyz
                                        * Referenced by: '<S246>/Iyz'
                                        */
  real_T Izz_Value;                    /* Expression: LD.Inertia.Izz
                                        * Referenced by: '<S246>/Izz'
                                        */
  real_T dIxxdt_Value;                 /* Expression: LD.Inertia.dIxx
                                        * Referenced by: '<S246>/dIxx//dt'
                                        */
  real_T dIxydt_Value;                 /* Expression: LD.Inertia.dIxy
                                        * Referenced by: '<S246>/dIxy//dt'
                                        */
  real_T dIxzdt_Value;                 /* Expression: LD.Inertia.dIxz
                                        * Referenced by: '<S246>/dIxz//dt'
                                        */
  real_T dIyydt_Value;                 /* Expression: LD.Inertia.dIyy
                                        * Referenced by: '<S246>/dIyy//dt'
                                        */
  real_T dIyzdt_Value;                 /* Expression: LD.Inertia.dIyz
                                        * Referenced by: '<S246>/dIyz//dt'
                                        */
  real_T dIzzdt_Value;                 /* Expression: LD.Inertia.dIzz
                                        * Referenced by: '<S246>/dIzz//dt'
                                        */
  real_T motor2position_Value_k[3];    /* Expression: LD.Propulsion.motor1.Position
                                        * Referenced by: '<S349>/motor2 position '
                                        */
  real_T Constant2_Value_pm;           /* Expression: 0
                                        * Referenced by: '<S349>/Constant2'
                                        */
  real_T Battery_Voltage_Value;        /* Expression: 29.6
                                        * Referenced by: '<S343>/Battery_Voltage'
                                        */
  real_T rendimiento_bateria_Value;    /* Expression: 0.98
                                        * Referenced by: '<S343>/rendimiento_bateria'
                                        */
  real_T eta_ESC_Value;                /* Expression: LD.Propulsion.ESC1.eta
                                        * Referenced by: '<S352>/eta_ESC'
                                        */
  real_T VoltagePrelookup_BreakpointsData[150];/* Expression: LD.Propulsion.motor1.Voltage
                                                * Referenced by: '<S353>/Voltage Prelookup'
                                                */
  real_T HeightPrelookup_BreakpointsData[5];/* Expression: LD.Propulsion.motor1.Height
                                             * Referenced by: '<S353>/Height Prelookup'
                                             */
  real_T VflightPrelookup_BreakpointsData_c[250];/* Expression: LD.Propulsion.motor1.Vflight
                                                  * Referenced by: '<S353>/Vflight Prelookup'
                                                  */
  real_T ThrustInterpolationUsingPrelookup_Table[187500];/* Expression: LD.Propulsion.motor1.Thrust
                                                          * Referenced by: '<S353>/Thrust Interpolation Using Prelookup'
                                                          */
  real_T Constant1_Value_fx;           /* Expression: 0
                                        * Referenced by: '<S349>/Constant1'
                                        */
  real_T RPMInterpolationUsingPrelookup_Table[187500];/* Expression: LD.Propulsion.motor1.RPM
                                                       * Referenced by: '<S353>/RPM Interpolation Using Prelookup'
                                                       */
  real_T Constant_Value_nn;            /* Expression: 1000
                                        * Referenced by: '<S354>/Constant'
                                        */
  real_T Constant1_Value_jd;           /* Expression: 0
                                        * Referenced by: '<S354>/Constant1'
                                        */
  real_T Switch_Threshold_c;           /* Expression: 0
                                        * Referenced by: '<S354>/Switch'
                                        */
  real_T Constant3_Value_nr;           /* Expression: 0
                                        * Referenced by: '<S349>/Constant3'
                                        */
  real_T Constant4_Value_ie;           /* Expression: 0
                                        * Referenced by: '<S349>/Constant4'
                                        */
  real_T motor2position_Value_i[3];    /* Expression: LD.Propulsion.motor2.Position
                                        * Referenced by: '<S363>/motor2 position '
                                        */
  real_T eta_ESC_Value_d;              /* Expression: LD.Propulsion.ESC2.eta
                                        * Referenced by: '<S366>/eta_ESC'
                                        */
  real_T VoltagePrelookup_BreakpointsData_k[150];/* Expression: LD.Propulsion.motor2.Voltage
                                                  * Referenced by: '<S367>/Voltage Prelookup'
                                                  */
  real_T HeightPrelookup_BreakpointsData_c[5];/* Expression: LD.Propulsion.motor2.Height
                                               * Referenced by: '<S367>/Height Prelookup'
                                               */
  real_T VflightPrelookup_BreakpointsData_p[250];/* Expression: LD.Propulsion.motor2.Vflight
                                                  * Referenced by: '<S367>/Vflight Prelookup'
                                                  */
  real_T ThrustInterpolationUsingPrelookup_Table_a[187500];/* Expression: LD.Propulsion.motor2.Thrust
                                                            * Referenced by: '<S367>/Thrust Interpolation Using Prelookup'
                                                            */
  real_T Constant2_Value_j;            /* Expression: 0
                                        * Referenced by: '<S363>/Constant2'
                                        */
  real_T Constant1_Value_c;            /* Expression: 0
                                        * Referenced by: '<S363>/Constant1'
                                        */
  real_T Constant4_Value_j;            /* Expression: 0
                                        * Referenced by: '<S363>/Constant4'
                                        */
  real_T Constant3_Value_j;            /* Expression: 0
                                        * Referenced by: '<S363>/Constant3'
                                        */
  real_T RPMInterpolationUsingPrelookup_Table_i[187500];/* Expression: LD.Propulsion.motor2.RPM
                                                         * Referenced by: '<S367>/RPM Interpolation Using Prelookup'
                                                         */
  real_T Constant_Value_b3;            /* Expression: 1000
                                        * Referenced by: '<S368>/Constant'
                                        */
  real_T Constant1_Value_lt;           /* Expression: 0
                                        * Referenced by: '<S368>/Constant1'
                                        */
  real_T Switch_Threshold_b;           /* Expression: 0
                                        * Referenced by: '<S368>/Switch'
                                        */
  real_T motor2position_Value_n[3];    /* Expression: LD.Propulsion.motor3.Position
                                        * Referenced by: '<S377>/motor2 position '
                                        */
  real_T eta_ESC_Value_n;              /* Expression: LD.Propulsion.ESC3.eta
                                        * Referenced by: '<S380>/eta_ESC'
                                        */
  real_T VoltagePrelookup_BreakpointsData_f[150];/* Expression: LD.Propulsion.motor3.Voltage
                                                  * Referenced by: '<S381>/Voltage Prelookup'
                                                  */
  real_T HeightPrelookup_BreakpointsData_b[5];/* Expression: LD.Propulsion.motor3.Height
                                               * Referenced by: '<S381>/Height Prelookup'
                                               */
  real_T VflightPrelookup_BreakpointsData_e[250];/* Expression: LD.Propulsion.motor3.Vflight
                                                  * Referenced by: '<S381>/Vflight Prelookup'
                                                  */
  real_T ThrustInterpolationUsingPrelookup_Table_k[187500];/* Expression: LD.Propulsion.motor3.Thrust
                                                            * Referenced by: '<S381>/Thrust Interpolation Using Prelookup'
                                                            */
  real_T Constant2_Value_d4;           /* Expression: 0
                                        * Referenced by: '<S377>/Constant2'
                                        */
  real_T Constant1_Value_cp;           /* Expression: 0
                                        * Referenced by: '<S377>/Constant1'
                                        */
  real_T Constant4_Value_o;            /* Expression: 0
                                        * Referenced by: '<S377>/Constant4'
                                        */
  real_T Constant3_Value_pv;           /* Expression: 0
                                        * Referenced by: '<S377>/Constant3'
                                        */
  real_T RPMInterpolationUsingPrelookup_Table_f[187500];/* Expression: LD.Propulsion.motor3.RPM
                                                         * Referenced by: '<S381>/RPM Interpolation Using Prelookup'
                                                         */
  real_T Constant_Value_hf;            /* Expression: 1000
                                        * Referenced by: '<S382>/Constant'
                                        */
  real_T Constant1_Value_ln;           /* Expression: 0
                                        * Referenced by: '<S382>/Constant1'
                                        */
  real_T Switch_Threshold_g;           /* Expression: 0
                                        * Referenced by: '<S382>/Switch'
                                        */
  real_T motor2position_Value_nw[3];   /* Expression: LD.Propulsion.motor4.Position
                                        * Referenced by: '<S391>/motor2 position '
                                        */
  real_T eta_ESC_Value_m;              /* Expression: LD.Propulsion.ESC4.eta
                                        * Referenced by: '<S394>/eta_ESC'
                                        */
  real_T VoltagePrelookup_BreakpointsData_m[150];/* Expression: LD.Propulsion.motor4.Voltage
                                                  * Referenced by: '<S395>/Voltage Prelookup'
                                                  */
  real_T HeightPrelookup_BreakpointsData_o[5];/* Expression: LD.Propulsion.motor4.Height
                                               * Referenced by: '<S395>/Height Prelookup'
                                               */
  real_T VflightPrelookup_BreakpointsData_jk[250];/* Expression: LD.Propulsion.motor4.Vflight
                                                   * Referenced by: '<S395>/Vflight Prelookup'
                                                   */
  real_T ThrustInterpolationUsingPrelookup_Table_i[187500];/* Expression: LD.Propulsion.motor4.Thrust
                                                            * Referenced by: '<S395>/Thrust Interpolation Using Prelookup'
                                                            */
  real_T Constant2_Value_l5;           /* Expression: 0
                                        * Referenced by: '<S391>/Constant2'
                                        */
  real_T Constant1_Value_ea;           /* Expression: 0
                                        * Referenced by: '<S391>/Constant1'
                                        */
  real_T Constant4_Value_dl;           /* Expression: 0
                                        * Referenced by: '<S391>/Constant4'
                                        */
  real_T Constant3_Value_nz;           /* Expression: 0
                                        * Referenced by: '<S391>/Constant3'
                                        */
  real_T RPMInterpolationUsingPrelookup_Table_p[187500];/* Expression: LD.Propulsion.motor4.RPM
                                                         * Referenced by: '<S395>/RPM Interpolation Using Prelookup'
                                                         */
  real_T Constant_Value_i2;            /* Expression: 1000
                                        * Referenced by: '<S396>/Constant'
                                        */
  real_T Constant1_Value_ct;           /* Expression: 0
                                        * Referenced by: '<S396>/Constant1'
                                        */
  real_T Switch_Threshold_a;           /* Expression: 0
                                        * Referenced by: '<S396>/Switch'
                                        */
  real_T motor2position_Value_l[3];    /* Expression: LD.Propulsion.motor5.Position
                                        * Referenced by: '<S405>/motor2 position '
                                        */
  real_T eta_ESC_Value_mv;             /* Expression: LD.Propulsion.ESC5.eta
                                        * Referenced by: '<S408>/eta_ESC'
                                        */
  real_T VoltagePrelookup_BreakpointsData_d[150];/* Expression: LD.Propulsion.motor5.Voltage
                                                  * Referenced by: '<S409>/Voltage Prelookup'
                                                  */
  real_T HeightPrelookup_BreakpointsData_h[5];/* Expression: LD.Propulsion.motor5.Height
                                               * Referenced by: '<S409>/Height Prelookup'
                                               */
  real_T VflightPrelookup_BreakpointsData_jw[250];/* Expression: LD.Propulsion.motor5.Vflight
                                                   * Referenced by: '<S409>/Vflight Prelookup'
                                                   */
  real_T ThrustInterpolationUsingPrelookup_Table_m[187500];/* Expression: LD.Propulsion.motor5.Thrust
                                                            * Referenced by: '<S409>/Thrust Interpolation Using Prelookup'
                                                            */
  real_T Constant2_Value_lh;           /* Expression: 0
                                        * Referenced by: '<S405>/Constant2'
                                        */
  real_T Constant1_Value_np;           /* Expression: 0
                                        * Referenced by: '<S405>/Constant1'
                                        */
  real_T Constant4_Value_kq;           /* Expression: 0
                                        * Referenced by: '<S405>/Constant4'
                                        */
  real_T Constant3_Value_ki;           /* Expression: 0
                                        * Referenced by: '<S405>/Constant3'
                                        */
  real_T RPMInterpolationUsingPrelookup_Table_k[187500];/* Expression: LD.Propulsion.motor5.RPM
                                                         * Referenced by: '<S409>/RPM Interpolation Using Prelookup'
                                                         */
  real_T Constant_Value_h4;            /* Expression: 1000
                                        * Referenced by: '<S410>/Constant'
                                        */
  real_T Constant1_Value_nj;           /* Expression: 0
                                        * Referenced by: '<S410>/Constant1'
                                        */
  real_T Switch_Threshold_g3;          /* Expression: 0
                                        * Referenced by: '<S410>/Switch'
                                        */
  real_T IC9_Value[3];                 /* Expression: [initialValues.Mx0 initialValues.My0 initialValues.Mz0]
                                        * Referenced by: '<S4>/IC9'
                                        */
  real_T dmassdtkgs_Value;             /* Expression: LD.Inertia.dmass
                                        * Referenced by: '<S246>/dmass//dt (kg//s)'
                                        */
  real_T Vre_xms_Value;                /* Expression: LD.Inertia.dmassdxb
                                        * Referenced by: '<S246>/Vre_x (m//s)'
                                        */
  real_T Vre_yms_Value;                /* Expression: LD.Inertia.dmassdyb
                                        * Referenced by: '<S246>/Vre_y (m//s)'
                                        */
  real_T Vre_zms_Value;                /* Expression: LD.Inertia.dmassdzb
                                        * Referenced by: '<S246>/Vre_z (m//s)'
                                        */
  real_T IC13_Value[3];                /* Expression: [-initialValues.D,0,-initialValues.L]
                                        * Referenced by: '<S4>/IC13'
                                        */
  real_T IC10_Value[3];                /* Expression: [initialValues.D,0,0]
                                        * Referenced by: '<S4>/IC10'
                                        */
  real_T IC6_Value[3];                 /* Expression: [initialValues.Fx0 initialValues.Fy0 initialValues.Fz0]
                                        * Referenced by: '<S4>/IC6'
                                        */
  real_T IC_Value_l[3];                /* Expression: [0 0 0]
                                        * Referenced by: '<S4>/IC'
                                        */
  real_T IC1_Value_n[3];               /* Expression: [0 0 0]
                                        * Referenced by: '<S4>/IC1'
                                        */
  real_T CurrentInterpolationUsingPrelookup_Table[187500];/* Expression: LD.Propulsion.motor1.Current
                                                           * Referenced by: '<S353>/Current Interpolation Using Prelookup'
                                                           */
  real_T CurrentInterpolationUsingPrelookup_Table_l[187500];/* Expression: LD.Propulsion.motor2.Current
                                                             * Referenced by: '<S367>/Current Interpolation Using Prelookup'
                                                             */
  real_T CurrentInterpolationUsingPrelookup_Table_e[187500];/* Expression: LD.Propulsion.motor3.Current
                                                             * Referenced by: '<S381>/Current Interpolation Using Prelookup'
                                                             */
  real_T CurrentInterpolationUsingPrelookup_Table_k[187500];/* Expression: LD.Propulsion.motor4.Current
                                                             * Referenced by: '<S395>/Current Interpolation Using Prelookup'
                                                             */
  real_T CurrentInterpolationUsingPrelookup_Table_lx[187500];/* Expression: LD.Propulsion.motor5.Current
                                                              * Referenced by: '<S409>/Current Interpolation Using Prelookup'
                                                              */
  real_T Constant_Value_is;            /* Expression: 1
                                        * Referenced by: '<S464>/Constant'
                                        */
  real_T UnitDelay_InitialCondition_b; /* Expression: 0
                                        * Referenced by: '<S464>/Unit Delay'
                                        */
  int32_T Constant1_Value_et;          /* Computed Parameter: Constant1_Value_et
                                        * Referenced by: '<S124>/Constant1'
                                        */
  int32_T Constant_Value_ho;           /* Computed Parameter: Constant_Value_ho
                                        * Referenced by: '<S125>/Constant'
                                        */
  int32_T Constant_Value_jj;           /* Computed Parameter: Constant_Value_jj
                                        * Referenced by: '<S123>/Constant'
                                        */
  int32_T Constant_Value_a;            /* Computed Parameter: Constant_Value_a
                                        * Referenced by: '<S134>/Constant'
                                        */
  int32_T Gain_Gain_ck;                /* Computed Parameter: Gain_Gain_ck
                                        * Referenced by: '<S134>/Gain'
                                        */
  int32_T Constant_Value_ci;           /* Computed Parameter: Constant_Value_ci
                                        * Referenced by: '<S137>/Constant'
                                        */
  int32_T Gain_Gain_ap;                /* Computed Parameter: Gain_Gain_ap
                                        * Referenced by: '<S136>/Gain'
                                        */
  int32_T Constant_Value_fe;           /* Computed Parameter: Constant_Value_fe
                                        * Referenced by: '<S140>/Constant'
                                        */
  int32_T Constant1_Value_at;          /* Computed Parameter: Constant1_Value_at
                                        * Referenced by: '<S140>/Constant1'
                                        */
  int32_T Constant1_Value_ge;          /* Computed Parameter: Constant1_Value_ge
                                        * Referenced by: '<S141>/Constant1'
                                        */
  int32_T Constant_Value_lu;           /* Computed Parameter: Constant_Value_lu
                                        * Referenced by: '<S139>/Constant'
                                        */
  int32_T Constant1_Value_b0;          /* Computed Parameter: Constant1_Value_b0
                                        * Referenced by: '<S138>/Constant1'
                                        */
  int32_T Gain_Gain_m;                 /* Computed Parameter: Gain_Gain_m
                                        * Referenced by: '<S138>/Gain'
                                        */
  int32_T Constant1_Value_eg;          /* Computed Parameter: Constant1_Value_eg
                                        * Referenced by: '<S142>/Constant1'
                                        */
  int32_T Constant_Value_a0;           /* Computed Parameter: Constant_Value_a0
                                        * Referenced by: '<S116>/Constant'
                                        */
  int32_T Constant_Value_g0;           /* Computed Parameter: Constant_Value_g0
                                        * Referenced by: '<S133>/Constant'
                                        */
  int32_T Gain_Gain_n;                 /* Computed Parameter: Gain_Gain_n
                                        * Referenced by: '<S133>/Gain'
                                        */
  int32_T Constant_Value_ha;           /* Computed Parameter: Constant_Value_ha
                                        * Referenced by: '<S143>/Constant'
                                        */
  int32_T Constant1_Value_e0;          /* Computed Parameter: Constant1_Value_e0
                                        * Referenced by: '<S143>/Constant1'
                                        */
  int32_T Constant_Value_et;           /* Computed Parameter: Constant_Value_et
                                        * Referenced by: '<S145>/Constant'
                                        */
  int32_T tc_old_Threshold;            /* Computed Parameter: tc_old_Threshold
                                        * Referenced by: '<S144>/tc_old'
                                        */
  int32_T Constant_Value_lrt;          /* Computed Parameter: Constant_Value_lrt
                                        * Referenced by: '<S115>/Constant'
                                        */
  int32_T Constant1_Value_dm;          /* Computed Parameter: Constant1_Value_dm
                                        * Referenced by: '<S115>/Constant1'
                                        */
  int32_T Constant_Value_fd;           /* Computed Parameter: Constant_Value_fd
                                        * Referenced by: '<S114>/Constant'
                                        */
  int32_T Constant_Value_dq;           /* Computed Parameter: Constant_Value_dq
                                        * Referenced by: '<S119>/Constant'
                                        */
  int32_T Gain_Gain_mq;                /* Computed Parameter: Gain_Gain_mq
                                        * Referenced by: '<S119>/Gain'
                                        */
  int32_T Constant_Value_bn;           /* Computed Parameter: Constant_Value_bn
                                        * Referenced by: '<S121>/Constant'
                                        */
  int32_T Constant_Value_io;           /* Computed Parameter: Constant_Value_io
                                        * Referenced by: '<S106>/Constant'
                                        */
  int32_T ForIterator_IterationLimit_f;/* Computed Parameter: ForIterator_IterationLimit_f
                                        * Referenced by: '<S106>/For Iterator'
                                        */
  int32_T arn_Threshold;               /* Computed Parameter: arn_Threshold
                                        * Referenced by: '<S106>/ar(n)'
                                        */
  int32_T offsetforupperindex1_Value;  /* Computed Parameter: offsetforupperindex1_Value
                                        * Referenced by: '<S340>/offset for upper index 1'
                                        */
  int32_T Constant_Value_g0l;          /* Computed Parameter: Constant_Value_g0l
                                        * Referenced by: '<S242>/Constant '
                                        */
  int32_T Constant14_Value_b;          /* Computed Parameter: Constant14_Value_b
                                        * Referenced by: '<S242>/Constant 14'
                                        */
  int32_T Constant2_Value_am;          /* Computed Parameter: Constant2_Value_am
                                        * Referenced by: '<S242>/Constant 2'
                                        */
  int32_T Constant4_Value_m;           /* Computed Parameter: Constant4_Value_m
                                        * Referenced by: '<S242>/Constant 4'
                                        */
  int32_T Constant6_Value_jg;          /* Computed Parameter: Constant6_Value_jg
                                        * Referenced by: '<S242>/Constant 6'
                                        */
  int32_T Constant8_Value_a;           /* Computed Parameter: Constant8_Value_a
                                        * Referenced by: '<S242>/Constant 8'
                                        */
  int32_T Constant10_Value_d;          /* Computed Parameter: Constant10_Value_d
                                        * Referenced by: '<S242>/Constant 10'
                                        */
  int32_T Constant12_Value_i;          /* Computed Parameter: Constant12_Value_i
                                        * Referenced by: '<S242>/Constant 12'
                                        */
  uint32_T TorqueInterpolationUsingPrelookup_maxIndex[2];/* Computed Parameter: TorqueInterpolationUsingPrelookup_maxIndex
                                                          * Referenced by: '<S354>/Torque Interpolation Using Prelookup'
                                                          */
  uint32_T TorqueInterpolationUsingPrelookup_maxIndex_e[2];/* Computed Parameter: TorqueInterpolationUsingPrelookup_maxIndex_e
                                                            * Referenced by: '<S368>/Torque Interpolation Using Prelookup'
                                                            */
  uint32_T TorqueInterpolationUsingPrelookup_maxIndex_b[2];/* Computed Parameter: TorqueInterpolationUsingPrelookup_maxIndex_b
                                                            * Referenced by: '<S382>/Torque Interpolation Using Prelookup'
                                                            */
  uint32_T TorqueInterpolationUsingPrelookup_maxIndex_g[2];/* Computed Parameter: TorqueInterpolationUsingPrelookup_maxIndex_g
                                                            * Referenced by: '<S396>/Torque Interpolation Using Prelookup'
                                                            */
  uint32_T TorqueInterpolationUsingPrelookup_maxIndex_c[2];/* Computed Parameter: TorqueInterpolationUsingPrelookup_maxIndex_c
                                                            * Referenced by: '<S410>/Torque Interpolation Using Prelookup'
                                                            */
  uint32_T MediumHighAltitudeIntensity_maxIndex[2];/* Computed Parameter: MediumHighAltitudeIntensity_maxIndex
                                                    * Referenced by: '<S73>/Medium//High Altitude Intensity'
                                                    */
  uint32_T CD0_dimSize[8];             /* Computed Parameter: CD0_dimSize
                                        * Referenced by: '<S253>/CD0'
                                        */
  uint32_T CD0_maxIndex[8];            /* Computed Parameter: CD0_maxIndex
                                        * Referenced by: '<S253>/CD0'
                                        */
  uint32_T CDalpha_dimSize[8];         /* Computed Parameter: CDalpha_dimSize
                                        * Referenced by: '<S253>/CDalpha'
                                        */
  uint32_T CDalpha_maxIndex[8];        /* Computed Parameter: CDalpha_maxIndex
                                        * Referenced by: '<S253>/CDalpha'
                                        */
  uint32_T CDalpha_dot_dimSize[8];     /* Computed Parameter: CDalpha_dot_dimSize
                                        * Referenced by: '<S253>/CDalpha_dot'
                                        */
  uint32_T CDalpha_dot_maxIndex[8];    /* Computed Parameter: CDalpha_dot_maxIndex
                                        * Referenced by: '<S253>/CDalpha_dot'
                                        */
  uint32_T CDq_dimSize[8];             /* Computed Parameter: CDq_dimSize
                                        * Referenced by: '<S253>/CDq'
                                        */
  uint32_T CDq_maxIndex[8];            /* Computed Parameter: CDq_maxIndex
                                        * Referenced by: '<S253>/CDq'
                                        */
  uint32_T CDdeltae_dimSize[8];        /* Computed Parameter: CDdeltae_dimSize
                                        * Referenced by: '<S253>/CDdeltae'
                                        */
  uint32_T CDdeltae_maxIndex[8];       /* Computed Parameter: CDdeltae_maxIndex
                                        * Referenced by: '<S253>/CDdeltae'
                                        */
  uint32_T CDdeltafr_dimSize[8];       /* Computed Parameter: CDdeltafr_dimSize
                                        * Referenced by: '<S253>/CDdeltafr'
                                        */
  uint32_T CDdeltafr_maxIndex[8];      /* Computed Parameter: CDdeltafr_maxIndex
                                        * Referenced by: '<S253>/CDdeltafr'
                                        */
  uint32_T CDdeltafrl_dimSize[8];      /* Computed Parameter: CDdeltafrl_dimSize
                                        * Referenced by: '<S253>/CDdeltafrl'
                                        */
  uint32_T CDdeltafrl_maxIndex[8];     /* Computed Parameter: CDdeltafrl_maxIndex
                                        * Referenced by: '<S253>/CDdeltafrl'
                                        */
  uint32_T CYbeta_dimSize[8];          /* Computed Parameter: CYbeta_dimSize
                                        * Referenced by: '<S255>/CYbeta'
                                        */
  uint32_T CYbeta_maxIndex[8];         /* Computed Parameter: CYbeta_maxIndex
                                        * Referenced by: '<S255>/CYbeta'
                                        */
  uint32_T CYp_dimSize[8];             /* Computed Parameter: CYp_dimSize
                                        * Referenced by: '<S255>/CYp'
                                        */
  uint32_T CYp_maxIndex[8];            /* Computed Parameter: CYp_maxIndex
                                        * Referenced by: '<S255>/CYp'
                                        */
  uint32_T CYr_dimSize[8];             /* Computed Parameter: CYr_dimSize
                                        * Referenced by: '<S255>/CYr'
                                        */
  uint32_T CYr_maxIndex[8];            /* Computed Parameter: CYr_maxIndex
                                        * Referenced by: '<S255>/CYr'
                                        */
  uint32_T CYdeltar_dimSize[8];        /* Computed Parameter: CYdeltar_dimSize
                                        * Referenced by: '<S255>/CYdeltar '
                                        */
  uint32_T CYdeltar_maxIndex[8];       /* Computed Parameter: CYdeltar_maxIndex
                                        * Referenced by: '<S255>/CYdeltar '
                                        */
  uint32_T CYdeltafr_dimSize[8];       /* Computed Parameter: CYdeltafr_dimSize
                                        * Referenced by: '<S255>/CYdeltafr'
                                        */
  uint32_T CYdeltafr_maxIndex[8];      /* Computed Parameter: CYdeltafr_maxIndex
                                        * Referenced by: '<S255>/CYdeltafr'
                                        */
  uint32_T CYdeltafl_dimSize[8];       /* Computed Parameter: CYdeltafl_dimSize
                                        * Referenced by: '<S255>/CYdeltafl'
                                        */
  uint32_T CYdeltafl_maxIndex[8];      /* Computed Parameter: CYdeltafl_maxIndex
                                        * Referenced by: '<S255>/CYdeltafl'
                                        */
  uint32_T CL0_dimSize[8];             /* Computed Parameter: CL0_dimSize
                                        * Referenced by: '<S254>/CL0'
                                        */
  uint32_T CLalpha_dimSize[8];         /* Computed Parameter: CLalpha_dimSize
                                        * Referenced by: '<S254>/CLalpha'
                                        */
  uint32_T CLalpha_maxIndex[8];        /* Computed Parameter: CLalpha_maxIndex
                                        * Referenced by: '<S254>/CLalpha'
                                        */
  uint32_T CLalpha_dot_dimSize[8];     /* Computed Parameter: CLalpha_dot_dimSize
                                        * Referenced by: '<S254>/CLalpha_dot'
                                        */
  uint32_T CLalpha_dot_maxIndex[8];    /* Computed Parameter: CLalpha_dot_maxIndex
                                        * Referenced by: '<S254>/CLalpha_dot'
                                        */
  uint32_T CLq_dimSize[8];             /* Computed Parameter: CLq_dimSize
                                        * Referenced by: '<S254>/CLq'
                                        */
  uint32_T CLq_maxIndex[8];            /* Computed Parameter: CLq_maxIndex
                                        * Referenced by: '<S254>/CLq'
                                        */
  uint32_T CLdeltae_dimSize[8];        /* Computed Parameter: CLdeltae_dimSize
                                        * Referenced by: '<S254>/CLdeltae'
                                        */
  uint32_T CLdeltae_maxIndex[8];       /* Computed Parameter: CLdeltae_maxIndex
                                        * Referenced by: '<S254>/CLdeltae'
                                        */
  uint32_T CLdeltafr_dimSize[8];       /* Computed Parameter: CLdeltafr_dimSize
                                        * Referenced by: '<S254>/CLdeltafr'
                                        */
  uint32_T CLdeltafr_maxIndex[8];      /* Computed Parameter: CLdeltafr_maxIndex
                                        * Referenced by: '<S254>/CLdeltafr'
                                        */
  uint32_T CLdeltafrl_dimSize[8];      /* Computed Parameter: CLdeltafrl_dimSize
                                        * Referenced by: '<S254>/CLdeltafrl'
                                        */
  uint32_T CLdeltafrl_maxIndex[8];     /* Computed Parameter: CLdeltafrl_maxIndex
                                        * Referenced by: '<S254>/CLdeltafrl'
                                        */
  uint32_T Clbeta_dimSize[8];          /* Computed Parameter: Clbeta_dimSize
                                        * Referenced by: '<S256>/Clbeta'
                                        */
  uint32_T Clbeta_maxIndex[8];         /* Computed Parameter: Clbeta_maxIndex
                                        * Referenced by: '<S256>/Clbeta'
                                        */
  uint32_T Clp_dimSize[8];             /* Computed Parameter: Clp_dimSize
                                        * Referenced by: '<S256>/Clp'
                                        */
  uint32_T Clp_maxIndex[8];            /* Computed Parameter: Clp_maxIndex
                                        * Referenced by: '<S256>/Clp'
                                        */
  uint32_T Clr_dimSize[8];             /* Computed Parameter: Clr_dimSize
                                        * Referenced by: '<S256>/Clr'
                                        */
  uint32_T Clr_maxIndex[8];            /* Computed Parameter: Clr_maxIndex
                                        * Referenced by: '<S256>/Clr'
                                        */
  uint32_T Cldeltar_dimSize[8];        /* Computed Parameter: Cldeltar_dimSize
                                        * Referenced by: '<S256>/Cldeltar '
                                        */
  uint32_T Cldeltar_maxIndex[8];       /* Computed Parameter: Cldeltar_maxIndex
                                        * Referenced by: '<S256>/Cldeltar '
                                        */
  uint32_T Cldeltafr_dimSize[8];       /* Computed Parameter: Cldeltafr_dimSize
                                        * Referenced by: '<S256>/Cldeltafr'
                                        */
  uint32_T Cldeltafr_maxIndex[8];      /* Computed Parameter: Cldeltafr_maxIndex
                                        * Referenced by: '<S256>/Cldeltafr'
                                        */
  uint32_T Cldeltafl_dimSize[8];       /* Computed Parameter: Cldeltafl_dimSize
                                        * Referenced by: '<S256>/Cldeltafl'
                                        */
  uint32_T Cldeltafl_maxIndex[8];      /* Computed Parameter: Cldeltafl_maxIndex
                                        * Referenced by: '<S256>/Cldeltafl'
                                        */
  uint32_T Cm0_dimSize[8];             /* Computed Parameter: Cm0_dimSize
                                        * Referenced by: '<S257>/Cm0'
                                        */
  uint32_T Cmalpha_dimSize[8];         /* Computed Parameter: Cmalpha_dimSize
                                        * Referenced by: '<S257>/Cmalpha'
                                        */
  uint32_T Cmalpha_maxIndex[8];        /* Computed Parameter: Cmalpha_maxIndex
                                        * Referenced by: '<S257>/Cmalpha'
                                        */
  uint32_T Cmalpha_dot_dimSize[8];     /* Computed Parameter: Cmalpha_dot_dimSize
                                        * Referenced by: '<S257>/Cmalpha_dot'
                                        */
  uint32_T Cmalpha_dot_maxIndex[8];    /* Computed Parameter: Cmalpha_dot_maxIndex
                                        * Referenced by: '<S257>/Cmalpha_dot'
                                        */
  uint32_T Cmq_dimSize[8];             /* Computed Parameter: Cmq_dimSize
                                        * Referenced by: '<S257>/Cmq'
                                        */
  uint32_T Cmq_maxIndex[8];            /* Computed Parameter: Cmq_maxIndex
                                        * Referenced by: '<S257>/Cmq'
                                        */
  uint32_T Cmdeltae_dimSize[8];        /* Computed Parameter: Cmdeltae_dimSize
                                        * Referenced by: '<S257>/Cmdeltae'
                                        */
  uint32_T Cmdeltae_maxIndex[8];       /* Computed Parameter: Cmdeltae_maxIndex
                                        * Referenced by: '<S257>/Cmdeltae'
                                        */
  uint32_T Cmdeltafr_dimSize[8];       /* Computed Parameter: Cmdeltafr_dimSize
                                        * Referenced by: '<S257>/Cmdeltafr'
                                        */
  uint32_T Cmdeltafr_maxIndex[8];      /* Computed Parameter: Cmdeltafr_maxIndex
                                        * Referenced by: '<S257>/Cmdeltafr'
                                        */
  uint32_T Cmdeltafl_dimSize[8];       /* Computed Parameter: Cmdeltafl_dimSize
                                        * Referenced by: '<S257>/Cmdeltafl'
                                        */
  uint32_T Cmdeltafl_maxIndex[8];      /* Computed Parameter: Cmdeltafl_maxIndex
                                        * Referenced by: '<S257>/Cmdeltafl'
                                        */
  uint32_T Cnbeta_dimSize[8];          /* Computed Parameter: Cnbeta_dimSize
                                        * Referenced by: '<S258>/Cnbeta'
                                        */
  uint32_T Cnbeta_maxIndex[8];         /* Computed Parameter: Cnbeta_maxIndex
                                        * Referenced by: '<S258>/Cnbeta'
                                        */
  uint32_T Cnp_dimSize[8];             /* Computed Parameter: Cnp_dimSize
                                        * Referenced by: '<S258>/Cnp'
                                        */
  uint32_T Cnp_maxIndex[8];            /* Computed Parameter: Cnp_maxIndex
                                        * Referenced by: '<S258>/Cnp'
                                        */
  uint32_T Cnr_dimSize[8];             /* Computed Parameter: Cnr_dimSize
                                        * Referenced by: '<S258>/Cnr'
                                        */
  uint32_T Cnr_maxIndex[8];            /* Computed Parameter: Cnr_maxIndex
                                        * Referenced by: '<S258>/Cnr'
                                        */
  uint32_T Cndeltar_dimSize[8];        /* Computed Parameter: Cndeltar_dimSize
                                        * Referenced by: '<S258>/Cndeltar '
                                        */
  uint32_T Cndeltar_maxIndex[8];       /* Computed Parameter: Cndeltar_maxIndex
                                        * Referenced by: '<S258>/Cndeltar '
                                        */
  uint32_T Cndeltafr_dimSize[8];       /* Computed Parameter: Cndeltafr_dimSize
                                        * Referenced by: '<S258>/Cndeltafr'
                                        */
  uint32_T Cndeltafr_maxIndex[8];      /* Computed Parameter: Cndeltafr_maxIndex
                                        * Referenced by: '<S258>/Cndeltafr'
                                        */
  uint32_T Cndeltafl_dimSize[8];       /* Computed Parameter: Cndeltafl_dimSize
                                        * Referenced by: '<S258>/Cndeltafl'
                                        */
  uint32_T Cndeltafl_maxIndex[8];      /* Computed Parameter: Cndeltafl_maxIndex
                                        * Referenced by: '<S258>/Cndeltafl'
                                        */
  uint32_T ThrustInterpolationUsingPrelookup_dimSize[3];/* Computed Parameter: ThrustInterpolationUsingPrelookup_dimSize
                                                         * Referenced by: '<S353>/Thrust Interpolation Using Prelookup'
                                                         */
  uint32_T ThrustInterpolationUsingPrelookup_maxIndex[3];/* Computed Parameter: ThrustInterpolationUsingPrelookup_maxIndex
                                                          * Referenced by: '<S353>/Thrust Interpolation Using Prelookup'
                                                          */
  uint32_T RPMInterpolationUsingPrelookup_dimSize[3];/* Computed Parameter: RPMInterpolationUsingPrelookup_dimSize
                                                      * Referenced by: '<S353>/RPM Interpolation Using Prelookup'
                                                      */
  uint32_T RPMInterpolationUsingPrelookup_maxIndex[3];/* Computed Parameter: RPMInterpolationUsingPrelookup_maxIndex
                                                       * Referenced by: '<S353>/RPM Interpolation Using Prelookup'
                                                       */
  uint32_T ThrustInterpolationUsingPrelookup_dimSize_c[3];/* Computed Parameter: ThrustInterpolationUsingPrelookup_dimSize_c
                                                           * Referenced by: '<S367>/Thrust Interpolation Using Prelookup'
                                                           */
  uint32_T ThrustInterpolationUsingPrelookup_maxIndex_g[3];/* Computed Parameter: ThrustInterpolationUsingPrelookup_maxIndex_g
                                                            * Referenced by: '<S367>/Thrust Interpolation Using Prelookup'
                                                            */
  uint32_T RPMInterpolationUsingPrelookup_dimSize_g[3];/* Computed Parameter: RPMInterpolationUsingPrelookup_dimSize_g
                                                        * Referenced by: '<S367>/RPM Interpolation Using Prelookup'
                                                        */
  uint32_T RPMInterpolationUsingPrelookup_maxIndex_m[3];/* Computed Parameter: RPMInterpolationUsingPrelookup_maxIndex_m
                                                         * Referenced by: '<S367>/RPM Interpolation Using Prelookup'
                                                         */
  uint32_T ThrustInterpolationUsingPrelookup_dimSize_b[3];/* Computed Parameter: ThrustInterpolationUsingPrelookup_dimSize_b
                                                           * Referenced by: '<S381>/Thrust Interpolation Using Prelookup'
                                                           */
  uint32_T ThrustInterpolationUsingPrelookup_maxIndex_o[3];/* Computed Parameter: ThrustInterpolationUsingPrelookup_maxIndex_o
                                                            * Referenced by: '<S381>/Thrust Interpolation Using Prelookup'
                                                            */
  uint32_T RPMInterpolationUsingPrelookup_dimSize_f[3];/* Computed Parameter: RPMInterpolationUsingPrelookup_dimSize_f
                                                        * Referenced by: '<S381>/RPM Interpolation Using Prelookup'
                                                        */
  uint32_T RPMInterpolationUsingPrelookup_maxIndex_mm[3];/* Computed Parameter: RPMInterpolationUsingPrelookup_maxIndex_mm
                                                          * Referenced by: '<S381>/RPM Interpolation Using Prelookup'
                                                          */
  uint32_T ThrustInterpolationUsingPrelookup_dimSize_bl[3];/* Computed Parameter: ThrustInterpolationUsingPrelookup_dimSize_bl
                                                            * Referenced by: '<S395>/Thrust Interpolation Using Prelookup'
                                                            */
  uint32_T ThrustInterpolationUsingPrelookup_maxIndex_a[3];/* Computed Parameter: ThrustInterpolationUsingPrelookup_maxIndex_a
                                                            * Referenced by: '<S395>/Thrust Interpolation Using Prelookup'
                                                            */
  uint32_T RPMInterpolationUsingPrelookup_dimSize_o[3];/* Computed Parameter: RPMInterpolationUsingPrelookup_dimSize_o
                                                        * Referenced by: '<S395>/RPM Interpolation Using Prelookup'
                                                        */
  uint32_T RPMInterpolationUsingPrelookup_maxIndex_c[3];/* Computed Parameter: RPMInterpolationUsingPrelookup_maxIndex_c
                                                         * Referenced by: '<S395>/RPM Interpolation Using Prelookup'
                                                         */
  uint32_T ThrustInterpolationUsingPrelookup_dimSize_n[3];/* Computed Parameter: ThrustInterpolationUsingPrelookup_dimSize_n
                                                           * Referenced by: '<S409>/Thrust Interpolation Using Prelookup'
                                                           */
  uint32_T ThrustInterpolationUsingPrelookup_maxIndex_p[3];/* Computed Parameter: ThrustInterpolationUsingPrelookup_maxIndex_p
                                                            * Referenced by: '<S409>/Thrust Interpolation Using Prelookup'
                                                            */
  uint32_T RPMInterpolationUsingPrelookup_dimSize_o5[3];/* Computed Parameter: RPMInterpolationUsingPrelookup_dimSize_o5
                                                         * Referenced by: '<S409>/RPM Interpolation Using Prelookup'
                                                         */
  uint32_T RPMInterpolationUsingPrelookup_maxIndex_h[3];/* Computed Parameter: RPMInterpolationUsingPrelookup_maxIndex_h
                                                         * Referenced by: '<S409>/RPM Interpolation Using Prelookup'
                                                         */
  uint32_T CurrentInterpolationUsingPrelookup_dimSize[3];/* Computed Parameter: CurrentInterpolationUsingPrelookup_dimSize
                                                          * Referenced by: '<S353>/Current Interpolation Using Prelookup'
                                                          */
  uint32_T CurrentInterpolationUsingPrelookup_maxIndex[3];/* Computed Parameter: CurrentInterpolationUsingPrelookup_maxIndex
                                                           * Referenced by: '<S353>/Current Interpolation Using Prelookup'
                                                           */
  uint32_T CurrentInterpolationUsingPrelookup_dimSize_e[3];/* Computed Parameter: CurrentInterpolationUsingPrelookup_dimSize_e
                                                            * Referenced by: '<S367>/Current Interpolation Using Prelookup'
                                                            */
  uint32_T CurrentInterpolationUsingPrelookup_maxIndex_h[3];/* Computed Parameter: CurrentInterpolationUsingPrelookup_maxIndex_h
                                                             * Referenced by: '<S367>/Current Interpolation Using Prelookup'
                                                             */
  uint32_T CurrentInterpolationUsingPrelookup_dimSize_p[3];/* Computed Parameter: CurrentInterpolationUsingPrelookup_dimSize_p
                                                            * Referenced by: '<S381>/Current Interpolation Using Prelookup'
                                                            */
  uint32_T CurrentInterpolationUsingPrelookup_maxIndex_o[3];/* Computed Parameter: CurrentInterpolationUsingPrelookup_maxIndex_o
                                                             * Referenced by: '<S381>/Current Interpolation Using Prelookup'
                                                             */
  uint32_T CurrentInterpolationUsingPrelookup_dimSize_k[3];/* Computed Parameter: CurrentInterpolationUsingPrelookup_dimSize_k
                                                            * Referenced by: '<S395>/Current Interpolation Using Prelookup'
                                                            */
  uint32_T CurrentInterpolationUsingPrelookup_maxIndex_c[3];/* Computed Parameter: CurrentInterpolationUsingPrelookup_maxIndex_c
                                                             * Referenced by: '<S395>/Current Interpolation Using Prelookup'
                                                             */
  uint32_T CurrentInterpolationUsingPrelookup_dimSize_i[3];/* Computed Parameter: CurrentInterpolationUsingPrelookup_dimSize_i
                                                            * Referenced by: '<S409>/Current Interpolation Using Prelookup'
                                                            */
  uint32_T CurrentInterpolationUsingPrelookup_maxIndex_a[3];/* Computed Parameter: CurrentInterpolationUsingPrelookup_maxIndex_a
                                                             * Referenced by: '<S409>/Current Interpolation Using Prelookup'
                                                             */
  uint8_T ManualSwitch7_CurrentSetting;/* Computed Parameter: ManualSwitch7_CurrentSetting
                                        * Referenced by: '<S164>/Manual Switch7'
                                        */
  uint8_T ManualSwitch_CurrentSetting; /* Computed Parameter: ManualSwitch_CurrentSetting
                                        * Referenced by: '<S164>/Manual Switch'
                                        */
  uint8_T ManualSwitch2_CurrentSetting;/* Computed Parameter: ManualSwitch2_CurrentSetting
                                        * Referenced by: '<S164>/Manual Switch2'
                                        */
  uint8_T ManualSwitch3_CurrentSetting;/* Computed Parameter: ManualSwitch3_CurrentSetting
                                        * Referenced by: '<S164>/Manual Switch3'
                                        */
  uint8_T ManualSwitch4_CurrentSetting;/* Computed Parameter: ManualSwitch4_CurrentSetting
                                        * Referenced by: '<S164>/Manual Switch4'
                                        */
  uint8_T ManualSwitch6_CurrentSetting;/* Computed Parameter: ManualSwitch6_CurrentSetting
                                        * Referenced by: '<S164>/Manual Switch6'
                                        */
  uint8_T ManualSwitch1_CurrentSetting;/* Computed Parameter: ManualSwitch1_CurrentSetting
                                        * Referenced by: '<S164>/Manual Switch1'
                                        */
  uint8_T ManualSwitch5_CurrentSetting;/* Computed Parameter: ManualSwitch5_CurrentSetting
                                        * Referenced by: '<S164>/Manual Switch5'
                                        */
  uint8_T ManualSwitch12_CurrentSetting;/* Computed Parameter: ManualSwitch12_CurrentSetting
                                         * Referenced by: '<S165>/Manual Switch12'
                                         */
  uint8_T ManualSwitch_CurrentSetting_f;/* Computed Parameter: ManualSwitch_CurrentSetting_f
                                         * Referenced by: '<S165>/Manual Switch'
                                         */
  uint8_T ManualSwitch1_CurrentSetting_h;/* Computed Parameter: ManualSwitch1_CurrentSetting_h
                                          * Referenced by: '<S165>/Manual Switch1'
                                          */
  uint8_T ManualSwitch6_CurrentSetting_n;/* Computed Parameter: ManualSwitch6_CurrentSetting_n
                                          * Referenced by: '<S165>/Manual Switch6'
                                          */
  uint8_T ManualSwitch5_CurrentSetting_d;/* Computed Parameter: ManualSwitch5_CurrentSetting_d
                                          * Referenced by: '<S165>/Manual Switch5'
                                          */
  uint8_T ManualSwitch3_CurrentSetting_i;/* Computed Parameter: ManualSwitch3_CurrentSetting_i
                                          * Referenced by: '<S165>/Manual Switch3'
                                          */
  uint8_T ManualSwitch11_CurrentSetting;/* Computed Parameter: ManualSwitch11_CurrentSetting
                                         * Referenced by: '<S166>/Manual Switch11'
                                         */
  uint8_T ManualSwitch10_CurrentSetting;/* Computed Parameter: ManualSwitch10_CurrentSetting
                                         * Referenced by: '<S166>/Manual Switch10'
                                         */
  uint8_T ManualSwitch2_CurrentSetting_c;/* Computed Parameter: ManualSwitch2_CurrentSetting_c
                                          * Referenced by: '<S166>/Manual Switch2'
                                          */
  uint8_T ManualSwitch4_CurrentSetting_b;/* Computed Parameter: ManualSwitch4_CurrentSetting_b
                                          * Referenced by: '<S166>/Manual Switch4'
                                          */
  uint8_T ManualSwitch13_CurrentSetting;/* Computed Parameter: ManualSwitch13_CurrentSetting
                                         * Referenced by: '<S166>/Manual Switch13'
                                         */
  uint8_T ManualSwitch7_CurrentSetting_f;/* Computed Parameter: ManualSwitch7_CurrentSetting_f
                                          * Referenced by: '<S166>/Manual Switch7'
                                          */
  uint8_T ManualSwitch8_CurrentSetting;/* Computed Parameter: ManualSwitch8_CurrentSetting
                                        * Referenced by: '<S166>/Manual Switch8'
                                        */
  uint8_T ManualSwitch12_CurrentSetting_k;/* Computed Parameter: ManualSwitch12_CurrentSetting_k
                                           * Referenced by: '<S166>/Manual Switch12'
                                           */
  uint8_T ManualSwitch_CurrentSetting_ff;/* Computed Parameter: ManualSwitch_CurrentSetting_ff
                                          * Referenced by: '<S166>/Manual Switch'
                                          */
  uint8_T ManualSwitch1_CurrentSetting_n;/* Computed Parameter: ManualSwitch1_CurrentSetting_n
                                          * Referenced by: '<S166>/Manual Switch1'
                                          */
  uint8_T ManualSwitch6_CurrentSetting_i;/* Computed Parameter: ManualSwitch6_CurrentSetting_i
                                          * Referenced by: '<S166>/Manual Switch6'
                                          */
  uint8_T ManualSwitch9_CurrentSetting;/* Computed Parameter: ManualSwitch9_CurrentSetting
                                        * Referenced by: '<S166>/Manual Switch9'
                                        */
  uint8_T ManualSwitch5_CurrentSetting_i;/* Computed Parameter: ManualSwitch5_CurrentSetting_i
                                          * Referenced by: '<S166>/Manual Switch5'
                                          */
  uint8_T ManualSwitch3_CurrentSetting_d;/* Computed Parameter: ManualSwitch3_CurrentSetting_d
                                          * Referenced by: '<S166>/Manual Switch3'
                                          */
  uint8_T ManualSwitch2_CurrentSetting_i;/* Computed Parameter: ManualSwitch2_CurrentSetting_i
                                          * Referenced by: '<S168>/Manual Switch2'
                                          */
  uint8_T ManualSwitch3_CurrentSetting_p;/* Computed Parameter: ManualSwitch3_CurrentSetting_p
                                          * Referenced by: '<S168>/Manual Switch3'
                                          */
  uint8_T ManualSwitch4_CurrentSetting_e;/* Computed Parameter: ManualSwitch4_CurrentSetting_e
                                          * Referenced by: '<S168>/Manual Switch4'
                                          */
  uint8_T ManualSwitch6_CurrentSetting_k;/* Computed Parameter: ManualSwitch6_CurrentSetting_k
                                          * Referenced by: '<S168>/Manual Switch6'
                                          */
  uint8_T ManualSwitch1_CurrentSetting_o;/* Computed Parameter: ManualSwitch1_CurrentSetting_o
                                          * Referenced by: '<S168>/Manual Switch1'
                                          */
  uint8_T ManualSwitch5_CurrentSetting_h;/* Computed Parameter: ManualSwitch5_CurrentSetting_h
                                          * Referenced by: '<S168>/Manual Switch5'
                                          */
  uint8_T ManualSwitch_CurrentSetting_o;/* Computed Parameter: ManualSwitch_CurrentSetting_o
                                         * Referenced by: '<S15>/Manual Switch'
                                         */
  uint8_T ManualSwitch1_CurrentSetting_g;/* Computed Parameter: ManualSwitch1_CurrentSetting_g
                                          * Referenced by: '<S15>/Manual Switch1'
                                          */
  uint8_T ManualSwitch_CurrentSetting_b;/* Computed Parameter: ManualSwitch_CurrentSetting_b
                                         * Referenced by: '<S5>/Manual Switch'
                                         */
  uint8_T ManualSwitch_CurrentSetting_oi;/* Computed Parameter: ManualSwitch_CurrentSetting_oi
                                          * Referenced by: '<S3>/Manual Switch'
                                          */
  P_JKFlipFlop_mainV03_56_T JKFlipFlop_jr;/* '<S476>/J-K Flip-Flop' */
  P_JKFlipFlop_mainV03_56_T JKFlipFlop_j;/* '<S475>/J-K Flip-Flop' */
  P_JKFlipFlop_mainV03_56_T JKFlipFlop_f;/* '<S474>/J-K Flip-Flop' */
  P_JKFlipFlop_mainV03_56_T JKFlipFlop_h;/* '<S473>/J-K Flip-Flop' */
  P_JKFlipFlop_mainV03_56_T JKFlipFlop_ej;/* '<S472>/J-K Flip-Flop' */
  P_JKFlipFlop_mainV03_56_T JKFlipFlop_n;/* '<S471>/J-K Flip-Flop' */
  P_JKFlipFlop_mainV03_56_T JKFlipFlop_b;/* '<S470>/J-K Flip-Flop' */
  P_JKFlipFlop_mainV03_56_T JKFlipFlop_p1;/* '<S469>/J-K Flip-Flop' */
  P_JKFlipFlop_mainV03_56_T JKFlipFlop_a;/* '<S468>/J-K Flip-Flop' */
  P_JKFlipFlop_mainV03_56_T JKFlipFlop_p;/* '<S467>/J-K Flip-Flop' */
  P_JKFlipFlop_mainV03_56_T JKFlipFlop_e;/* '<S466>/J-K Flip-Flop' */
  P_JKFlipFlop_mainV03_56_T JKFlipFlop;/* '<S465>/J-K Flip-Flop' */
  P_Distanceintogusty_mainV03_56_T Distanceintogustz;/* '<S47>/Distance into gust (z)' */
  P_Distanceintogusty_mainV03_56_T Distanceintogusty;/* '<S47>/Distance into gust (y)' */
};

/* Real-time Model Data Structure */
struct tag_RTM_mainV03_56_T {
  const char_T *path;
  const char_T *modelName;
  struct SimStruct_tag * *childSfunctions;
  const char_T *errorStatus;
  SS_SimMode simMode;
  RTWLogInfo *rtwLogInfo;
  RTWExtModeInfo *extModeInfo;
  RTWSolverInfo solverInfo;
  RTWSolverInfo *solverInfoPtr;
  void *sfcnInfo;

  /*
   * NonInlinedSFcns:
   * The following substructure contains information regarding
   * non-inlined s-functions used in the model.
   */
  struct {
    RTWSfcnInfo sfcnInfo;
    time_T *taskTimePtrs[3];
    SimStruct childSFunctions[1];
    SimStruct *childSFunctionPtrs[1];
    struct _ssBlkInfo2 blkInfo2[1];
    struct _ssSFcnModelMethods2 methods2[1];
    struct _ssSFcnModelMethods3 methods3[1];
    struct _ssSFcnModelMethods4 methods4[1];
    struct _ssStatesInfo2 statesInfo2[1];
    ssPeriodicStatesInfo periodicStatesInfo[1];
    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortOutputs outputPortInfo[3];
      uint_T attribs[3];
      mxArray *params[3];
      struct _ssDWorkRecord dWork[2];
      struct _ssDWorkAuxRecord dWorkAux[2];
    } Sfcn0;
  } NonInlinedSFcns;

  void *blockIO;
  const void *constBlockIO;
  void *defaultParam;
  ZCSigState *prevZCSigState;
  real_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  void *zcSignalValues;
  void *inputs;
  void *outputs;
  boolean_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T blkStateChange;
  real_T odeX0[140];
  real_T odeF0[140];
  real_T odeX1START[140];
  real_T odeF1[140];
  real_T odeDELTA[140];
  real_T odeE[4*140];
  real_T odeFAC[140];
  real_T odeDFDX[140*140];
  real_T odeW[140*140];
  int_T odePIVOTS[140];
  real_T odeXTMP[140];
  real_T odeZTMP[140];
  ODE14X_IntgData intgData;
  void *dwork;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    uint32_T checksums[4];
    uint32_T options;
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numU;
    int_T numY;
    int_T numSampTimes;
    int_T numBlocks;
    int_T numBlockIO;
    int_T numBlockPrms;
    int_T numDwork;
    int_T numSFcnPrms;
    int_T numSFcns;
    int_T numIports;
    int_T numOports;
    int_T numNonSampZCs;
    int_T sysDirFeedThru;
    int_T rtwGenSfcn;
  } Sizes;

  /*
   * SpecialInfo:
   * The following substructure contains special information
   * related to other components that are dependent on RTW.
   */
  struct {
    const void *mappingInfo;
    void *xpcData;
  } SpecialInfo;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    time_T stepSize;
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    time_T stepSize1;
    uint32_T clockTick2;
    uint32_T clockTickH2;
    time_T stepSize2;
    struct {
      uint8_T TID[3];
    } TaskCounters;

    time_T tStart;
    time_T tFinal;
    time_T timeOfLastOutput;
    void *timingData;
    real_T *varNextHitTimesList;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *sampleTimes;
    time_T *offsetTimes;
    int_T *sampleTimeTaskIDPtr;
    int_T *sampleHits;
    int_T *perTaskSampleHits;
    time_T *t;
    time_T sampleTimesArray[3];
    time_T offsetTimesArray[3];
    int_T sampleTimeTaskIDArray[3];
    int_T sampleHitArray[3];
    int_T perTaskSampleHitsArray[9];
    time_T tArray[3];
  } Timing;
};

/* Block parameters (auto storage) */
extern P_mainV03_56_T mainV03_56_P;

/* Block signals (auto storage) */
extern B_mainV03_56_T mainV03_56_B;

/* Continuous states (auto storage) */
extern X_mainV03_56_T mainV03_56_X;

/* Block states (auto storage) */
extern DW_mainV03_56_T mainV03_56_DW;

/* Model entry point functions */
extern void mainV03_56_initialize(void);
extern void mainV03_56_output(void);
extern void mainV03_56_update(void);
extern void mainV03_56_terminate(void);

/*====================*
 * External functions *
 *====================*/
extern mainV03_56_rtModel *mainV03_56(void);
extern void MdlInitializeSizes(void);
extern void MdlInitializeSampleTimes(void);
extern void MdlInitialize(void);
extern void MdlStart(void);
extern void MdlOutputs(int_T tid);
extern void MdlUpdate(int_T tid);
extern void MdlTerminate(void);

/* Real-time Model object */
extern RT_MODEL_mainV03_56_T *const mainV03_56_M;

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'mainV03_56'
 * '<S1>'   : 'mainV03_56/Actuators'
 * '<S2>'   : 'mainV03_56/Environment'
 * '<S3>'   : 'mainV03_56/FCS'
 * '<S4>'   : 'mainV03_56/LIBIS Plant'
 * '<S5>'   : 'mainV03_56/Pilot Commands'
 * '<S6>'   : 'mainV03_56/Sensors'
 * '<S7>'   : 'mainV03_56/Visualization'
 * '<S8>'   : 'mainV03_56/Actuators/Feedthrough'
 * '<S9>'   : 'mainV03_56/Actuators/Linear 1st Order'
 * '<S10>'  : 'mainV03_56/Actuators/Linear 2nd Order'
 * '<S11>'  : 'mainV03_56/Actuators/Nonlinear 2nd Order'
 * '<S12>'  : 'mainV03_56/Environment/COESA Atmosphere Model'
 * '<S13>'  : 'mainV03_56/Environment/Terrain'
 * '<S14>'  : 'mainV03_56/Environment/WGS84 Gravity Model '
 * '<S15>'  : 'mainV03_56/Environment/Wind Models'
 * '<S16>'  : 'mainV03_56/Environment/World Magnetic Model 2015'
 * '<S17>'  : 'mainV03_56/Environment/COESA Atmosphere Model/Density Conversion'
 * '<S18>'  : 'mainV03_56/Environment/COESA Atmosphere Model/Length Conversion'
 * '<S19>'  : 'mainV03_56/Environment/COESA Atmosphere Model/Pressure Conversion'
 * '<S20>'  : 'mainV03_56/Environment/COESA Atmosphere Model/Temperature Conversion'
 * '<S21>'  : 'mainV03_56/Environment/COESA Atmosphere Model/Velocity Conversion'
 * '<S22>'  : 'mainV03_56/Environment/Terrain/3x3 Cross Product'
 * '<S23>'  : 'mainV03_56/Environment/Terrain/3x3 Cross Product1'
 * '<S24>'  : 'mainV03_56/Environment/Terrain/3x3 Cross Product2'
 * '<S25>'  : 'mainV03_56/Environment/Terrain/3x3 Cross Product3'
 * '<S26>'  : 'mainV03_56/Environment/Terrain/3x3 Cross Product4'
 * '<S27>'  : 'mainV03_56/Environment/Terrain/3x3 Cross Product5'
 * '<S28>'  : 'mainV03_56/Environment/Terrain/Reaction'
 * '<S29>'  : 'mainV03_56/Environment/Terrain/Reaction1'
 * '<S30>'  : 'mainV03_56/Environment/Terrain/Reaction2'
 * '<S31>'  : 'mainV03_56/Environment/Terrain/3x3 Cross Product/Subsystem'
 * '<S32>'  : 'mainV03_56/Environment/Terrain/3x3 Cross Product/Subsystem1'
 * '<S33>'  : 'mainV03_56/Environment/Terrain/3x3 Cross Product1/Subsystem'
 * '<S34>'  : 'mainV03_56/Environment/Terrain/3x3 Cross Product1/Subsystem1'
 * '<S35>'  : 'mainV03_56/Environment/Terrain/3x3 Cross Product2/Subsystem'
 * '<S36>'  : 'mainV03_56/Environment/Terrain/3x3 Cross Product2/Subsystem1'
 * '<S37>'  : 'mainV03_56/Environment/Terrain/3x3 Cross Product3/Subsystem'
 * '<S38>'  : 'mainV03_56/Environment/Terrain/3x3 Cross Product3/Subsystem1'
 * '<S39>'  : 'mainV03_56/Environment/Terrain/3x3 Cross Product4/Subsystem'
 * '<S40>'  : 'mainV03_56/Environment/Terrain/3x3 Cross Product4/Subsystem1'
 * '<S41>'  : 'mainV03_56/Environment/Terrain/3x3 Cross Product5/Subsystem'
 * '<S42>'  : 'mainV03_56/Environment/Terrain/3x3 Cross Product5/Subsystem1'
 * '<S43>'  : 'mainV03_56/Environment/WGS84 Gravity Model /Acceleration Conversion'
 * '<S44>'  : 'mainV03_56/Environment/WGS84 Gravity Model /Angle Conversion'
 * '<S45>'  : 'mainV03_56/Environment/WGS84 Gravity Model /Length Conversion'
 * '<S46>'  : 'mainV03_56/Environment/WGS84 Gravity Model /Velocity Conversion2'
 * '<S47>'  : 'mainV03_56/Environment/Wind Models/Discrete Wind Gust Model'
 * '<S48>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))'
 * '<S49>'  : 'mainV03_56/Environment/Wind Models/Wind Shear Model'
 * '<S50>'  : 'mainV03_56/Environment/Wind Models/Discrete Wind Gust Model/Distance into gust (x)'
 * '<S51>'  : 'mainV03_56/Environment/Wind Models/Discrete Wind Gust Model/Distance into gust (y)'
 * '<S52>'  : 'mainV03_56/Environment/Wind Models/Discrete Wind Gust Model/Distance into gust (z)'
 * '<S53>'  : 'mainV03_56/Environment/Wind Models/Discrete Wind Gust Model/Velocity Conversion'
 * '<S54>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Angle Conversion'
 * '<S55>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Filters on angular rates'
 * '<S56>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Filters on velocities'
 * '<S57>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Length Conversion'
 * '<S58>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Length Conversion1'
 * '<S59>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/RMS turbulence  intensities'
 * '<S60>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Select angular rates'
 * '<S61>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Select velocities'
 * '<S62>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Turbulence scale lengths'
 * '<S63>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Velocity Conversion'
 * '<S64>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Velocity Conversion2'
 * '<S65>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/White Noise (p)'
 * '<S66>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/White Noise (u,v,w)'
 * '<S67>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Filters on angular rates/Hpgw'
 * '<S68>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Filters on angular rates/Hqgw'
 * '<S69>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Filters on angular rates/Hrgw'
 * '<S70>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Filters on velocities/Hugw(z)'
 * '<S71>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Filters on velocities/Hvgw(z)'
 * '<S72>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Filters on velocities/Hwgw(z)'
 * '<S73>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/RMS turbulence  intensities/High Altitude Intensity'
 * '<S74>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/RMS turbulence  intensities/Low Altitude Intensity'
 * '<S75>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Select angular rates/Interpolate  rates'
 * '<S76>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Select angular rates/Low altitude  rates'
 * '<S77>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Select angular rates/Medium//High  altitude rates'
 * '<S78>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Select angular rates/Merge Subsystems'
 * '<S79>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Select angular rates/Interpolate  rates/wind to body transformation'
 * '<S80>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Select angular rates/Interpolate  rates/wind to body transformation/convert to earth coords'
 * '<S81>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Select angular rates/Low altitude  rates/wind to body transformation'
 * '<S82>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Select angular rates/Low altitude  rates/wind to body transformation/convert to earth coords'
 * '<S83>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Select velocities/Interpolate  velocities'
 * '<S84>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Select velocities/Low altitude  velocities'
 * '<S85>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Select velocities/Medium//High  altitude velocities'
 * '<S86>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Select velocities/Merge Subsystems'
 * '<S87>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Select velocities/Interpolate  velocities/wind to body transformation'
 * '<S88>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Select velocities/Interpolate  velocities/wind to body transformation/convert to earth coords'
 * '<S89>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Select velocities/Low altitude  velocities/wind to body transformation'
 * '<S90>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Select velocities/Low altitude  velocities/wind to body transformation/convert to earth coords'
 * '<S91>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Turbulence scale lengths/Low altitude scale length'
 * '<S92>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Turbulence scale lengths/Medium//High altitude scale length'
 * '<S93>'  : 'mainV03_56/Environment/Wind Models/Dryden Wind Turbulence Model  (Discrete (+q +r))/Turbulence scale lengths/Medium//High altitude scale length/Length Conversion'
 * '<S94>'  : 'mainV03_56/Environment/Wind Models/Wind Shear Model/Angle Conversion'
 * '<S95>'  : 'mainV03_56/Environment/Wind Models/Wind Shear Model/Length Conversion'
 * '<S96>'  : 'mainV03_56/Environment/World Magnetic Model 2015/Check Altitude'
 * '<S97>'  : 'mainV03_56/Environment/World Magnetic Model 2015/Check Latitude'
 * '<S98>'  : 'mainV03_56/Environment/World Magnetic Model 2015/Check Longitude'
 * '<S99>'  : 'mainV03_56/Environment/World Magnetic Model 2015/Compute x,y,z, and h components of magnetic field'
 * '<S100>' : 'mainV03_56/Environment/World Magnetic Model 2015/Is time within model limits'
 * '<S101>' : 'mainV03_56/Environment/World Magnetic Model 2015/Length Conversion'
 * '<S102>' : 'mainV03_56/Environment/World Magnetic Model 2015/MagField Conversion'
 * '<S103>' : 'mainV03_56/Environment/World Magnetic Model 2015/MagField Conversion1'
 * '<S104>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag'
 * '<S105>' : 'mainV03_56/Environment/World Magnetic Model 2015/Compute x,y,z, and h components of magnetic field/Angle Conversion'
 * '<S106>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates'
 * '<S107>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Convert from geodetic to  spherical coordinates'
 * '<S108>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Convert from geodetic to  spherical coordinates '
 * '<S109>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Get Cosine and Sine  of Latitude and Longitude'
 * '<S110>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Has altitude or latitude changed'
 * '<S111>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Has longitude changed '
 * '<S112>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Has time changed'
 * '<S113>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Rotate magnetic vector components  to geodetic from spherical and  compute declination, inclination  and total intensity'
 * '<S114>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem'
 * '<S115>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Accumulate terms of the  spherical harmonic expansion'
 * '<S116>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations'
 * '<S117>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Time adjust the gauss coefficients'
 * '<S118>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Accumulate terms of the  spherical harmonic expansion/Special case - North//South Geographic Pole'
 * '<S119>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Accumulate terms of the  spherical harmonic expansion/calculate  index'
 * '<S120>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Accumulate terms of the  spherical harmonic expansion/calculate temp values'
 * '<S121>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Accumulate terms of the  spherical harmonic expansion/special case'
 * '<S122>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Accumulate terms of the  spherical harmonic expansion/Special case - North//South Geographic Pole/If Action Subsystem1'
 * '<S123>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Accumulate terms of the  spherical harmonic expansion/Special case - North//South Geographic Pole/If Action Subsystem2'
 * '<S124>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Accumulate terms of the  spherical harmonic expansion/Special case - North//South Geographic Pole/If Action Subsystem2/calculate  indices'
 * '<S125>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Accumulate terms of the  spherical harmonic expansion/Special case - North//South Geographic Pole/If Action Subsystem2/calculate  row and column'
 * '<S126>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Accumulate terms of the  spherical harmonic expansion/calculate temp values/If Action Subsystem'
 * '<S127>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Accumulate terms of the  spherical harmonic expansion/calculate temp values/If Action Subsystem1'
 * '<S128>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Accumulate terms of the  spherical harmonic expansion/calculate temp values/If Action Subsystem1/m,n'
 * '<S129>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Accumulate terms of the  spherical harmonic expansion/calculate temp values/If Action Subsystem1/n,m-1'
 * '<S130>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations/If Action Subsystem'
 * '<S131>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations/If Action Subsystem1'
 * '<S132>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations/If Action Subsystem2'
 * '<S133>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations/calculate  index'
 * '<S134>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations/If Action Subsystem/calculate  index'
 * '<S135>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations/If Action Subsystem/calculate  row and column'
 * '<S136>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations/If Action Subsystem1/calculate  index'
 * '<S137>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations/If Action Subsystem1/calculate  row and column'
 * '<S138>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations/If Action Subsystem2/calculate  indices'
 * '<S139>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations/If Action Subsystem2/calculate  row and column1'
 * '<S140>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations/If Action Subsystem2/calculate  row and column2'
 * '<S141>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations/If Action Subsystem2/m<n-2'
 * '<S142>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations/If Action Subsystem2/m<n-2 '
 * '<S143>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Time adjust the gauss coefficients/If Action Subsystem'
 * '<S144>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Time adjust the gauss coefficients/if (m~=0)'
 * '<S145>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Time adjust the gauss coefficients/if (m~=0)/If Action Subsystem1'
 * '<S146>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Time adjust the gauss coefficients/if (m~=0)/If Action Subsystem2'
 * '<S147>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Convert from geodetic to  spherical coordinates/calculate ca'
 * '<S148>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Convert from geodetic to  spherical coordinates/calculate ct'
 * '<S149>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Convert from geodetic to  spherical coordinates/calculate d'
 * '<S150>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Convert from geodetic to  spherical coordinates/calculate q'
 * '<S151>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Convert from geodetic to  spherical coordinates/calculate q2'
 * '<S152>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Convert from geodetic to  spherical coordinates/calculate r2'
 * '<S153>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Convert from geodetic to  spherical coordinates/calculate sa'
 * '<S154>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Convert from geodetic to  spherical coordinates/calculate st'
 * '<S155>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Convert from geodetic to  spherical coordinates /For Iterator Subsystem'
 * '<S156>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Get Cosine and Sine  of Latitude and Longitude/Angle Conversion2'
 * '<S157>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Rotate magnetic vector components  to geodetic from spherical and  compute declination, inclination  and total intensity/Calculate bx'
 * '<S158>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Rotate magnetic vector components  to geodetic from spherical and  compute declination, inclination  and total intensity/Calculate by'
 * '<S159>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Rotate magnetic vector components  to geodetic from spherical and  compute declination, inclination  and total intensity/Calculate bz'
 * '<S160>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Rotate magnetic vector components  to geodetic from spherical and  compute declination, inclination  and total intensity/Compute declination, inclination,  and total intensity'
 * '<S161>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Rotate magnetic vector components  to geodetic from spherical and  compute declination, inclination  and total intensity/Compute declination, inclination,  and total intensity/Angle Conversion'
 * '<S162>' : 'mainV03_56/Environment/World Magnetic Model 2015/geomag/Rotate magnetic vector components  to geodetic from spherical and  compute declination, inclination  and total intensity/Compute declination, inclination,  and total intensity/Angle Conversion1'
 * '<S163>' : 'mainV03_56/FCS/Fixed-Wing - Cruise'
 * '<S164>' : 'mainV03_56/FCS/Fixed-Wing Climb'
 * '<S165>' : 'mainV03_56/FCS/Quadcopter'
 * '<S166>' : 'mainV03_56/FCS/Quadcopter --> Fixed-Wing'
 * '<S167>' : 'mainV03_56/FCS/Select Active FCS'
 * '<S168>' : 'mainV03_56/FCS/TakeOff'
 * '<S169>' : 'mainV03_56/FCS/Fixed-Wing Climb/P Height'
 * '<S170>' : 'mainV03_56/FCS/Fixed-Wing Climb/PID Controller'
 * '<S171>' : 'mainV03_56/FCS/Fixed-Wing Climb/PID Controller1'
 * '<S172>' : 'mainV03_56/FCS/Fixed-Wing Climb/PID Controller2'
 * '<S173>' : 'mainV03_56/FCS/Fixed-Wing Climb/PID Controller3'
 * '<S174>' : 'mainV03_56/FCS/Fixed-Wing Climb/PID Controller4'
 * '<S175>' : 'mainV03_56/FCS/Fixed-Wing Climb/PID Controller5'
 * '<S176>' : 'mainV03_56/FCS/Fixed-Wing Climb/PID Controller7'
 * '<S177>' : 'mainV03_56/FCS/Fixed-Wing Climb/Rotation Angles to Direction Cosine Matrix'
 * '<S178>' : 'mainV03_56/FCS/Fixed-Wing Climb/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S179>' : 'mainV03_56/FCS/Quadcopter/ESTE2'
 * '<S180>' : 'mainV03_56/FCS/Quadcopter/PID Altitude Acceleration'
 * '<S181>' : 'mainV03_56/FCS/Quadcopter/PID Altitude Rate'
 * '<S182>' : 'mainV03_56/FCS/Quadcopter/PID Controller'
 * '<S183>' : 'mainV03_56/FCS/Quadcopter/PID Controller1'
 * '<S184>' : 'mainV03_56/FCS/Quadcopter/PID Height Controller1'
 * '<S185>' : 'mainV03_56/FCS/Quadcopter/PID Lateral Position'
 * '<S186>' : 'mainV03_56/FCS/Quadcopter/PID Pitch Rate2'
 * '<S187>' : 'mainV03_56/FCS/Quadcopter/PID Roll Rate'
 * '<S188>' : 'mainV03_56/FCS/Quadcopter/PID Yaw Rate1'
 * '<S189>' : 'mainV03_56/FCS/Quadcopter/Proportional Altitude'
 * '<S190>' : 'mainV03_56/FCS/Quadcopter/Proportional Pitch2'
 * '<S191>' : 'mainV03_56/FCS/Quadcopter/Proportional Roll'
 * '<S192>' : 'mainV03_56/FCS/Quadcopter/Proportional Yaw1'
 * '<S193>' : 'mainV03_56/FCS/Quadcopter/Quad actuator commands'
 * '<S194>' : 'mainV03_56/FCS/Quadcopter/Rotation Angles to Direction Cosine Matrix'
 * '<S195>' : 'mainV03_56/FCS/Quadcopter/Quad actuator commands/Altitude'
 * '<S196>' : 'mainV03_56/FCS/Quadcopter/Quad actuator commands/Pitch'
 * '<S197>' : 'mainV03_56/FCS/Quadcopter/Quad actuator commands/Roll'
 * '<S198>' : 'mainV03_56/FCS/Quadcopter/Quad actuator commands/Yaw'
 * '<S199>' : 'mainV03_56/FCS/Quadcopter/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S200>' : 'mainV03_56/FCS/Quadcopter --> Fixed-Wing/ESTE2'
 * '<S201>' : 'mainV03_56/FCS/Quadcopter --> Fixed-Wing/P Height1'
 * '<S202>' : 'mainV03_56/FCS/Quadcopter --> Fixed-Wing/PID Altitude Acceleration'
 * '<S203>' : 'mainV03_56/FCS/Quadcopter --> Fixed-Wing/PID Altitude Rate'
 * '<S204>' : 'mainV03_56/FCS/Quadcopter --> Fixed-Wing/PID Controller'
 * '<S205>' : 'mainV03_56/FCS/Quadcopter --> Fixed-Wing/PID Controller1'
 * '<S206>' : 'mainV03_56/FCS/Quadcopter --> Fixed-Wing/PID Controller4'
 * '<S207>' : 'mainV03_56/FCS/Quadcopter --> Fixed-Wing/PID Controller5'
 * '<S208>' : 'mainV03_56/FCS/Quadcopter --> Fixed-Wing/PID Controller6'
 * '<S209>' : 'mainV03_56/FCS/Quadcopter --> Fixed-Wing/PID Controller7'
 * '<S210>' : 'mainV03_56/FCS/Quadcopter --> Fixed-Wing/PID Controller8'
 * '<S211>' : 'mainV03_56/FCS/Quadcopter --> Fixed-Wing/PID Controller9'
 * '<S212>' : 'mainV03_56/FCS/Quadcopter --> Fixed-Wing/PID Height Controller1'
 * '<S213>' : 'mainV03_56/FCS/Quadcopter --> Fixed-Wing/PID Lateral Position'
 * '<S214>' : 'mainV03_56/FCS/Quadcopter --> Fixed-Wing/PID Pitch Rate2'
 * '<S215>' : 'mainV03_56/FCS/Quadcopter --> Fixed-Wing/PID Roll Rate'
 * '<S216>' : 'mainV03_56/FCS/Quadcopter --> Fixed-Wing/PID Yaw Rate1'
 * '<S217>' : 'mainV03_56/FCS/Quadcopter --> Fixed-Wing/Proportional Altitude'
 * '<S218>' : 'mainV03_56/FCS/Quadcopter --> Fixed-Wing/Proportional Pitch2'
 * '<S219>' : 'mainV03_56/FCS/Quadcopter --> Fixed-Wing/Proportional Roll'
 * '<S220>' : 'mainV03_56/FCS/Quadcopter --> Fixed-Wing/Proportional Yaw1'
 * '<S221>' : 'mainV03_56/FCS/Quadcopter --> Fixed-Wing/Quad actuator commands'
 * '<S222>' : 'mainV03_56/FCS/Quadcopter --> Fixed-Wing/Rotation Angles to Direction Cosine Matrix'
 * '<S223>' : 'mainV03_56/FCS/Quadcopter --> Fixed-Wing/Quad actuator commands/Altitude'
 * '<S224>' : 'mainV03_56/FCS/Quadcopter --> Fixed-Wing/Quad actuator commands/Pitch'
 * '<S225>' : 'mainV03_56/FCS/Quadcopter --> Fixed-Wing/Quad actuator commands/Roll'
 * '<S226>' : 'mainV03_56/FCS/Quadcopter --> Fixed-Wing/Quad actuator commands/Yaw'
 * '<S227>' : 'mainV03_56/FCS/Quadcopter --> Fixed-Wing/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S228>' : 'mainV03_56/FCS/Select Active FCS/Compare To Constant'
 * '<S229>' : 'mainV03_56/FCS/Select Active FCS/Compare To Constant1'
 * '<S230>' : 'mainV03_56/FCS/Select Active FCS/Compare To Constant2'
 * '<S231>' : 'mainV03_56/FCS/Select Active FCS/Compare To Constant3'
 * '<S232>' : 'mainV03_56/FCS/Select Active FCS/Compare To Constant4'
 * '<S233>' : 'mainV03_56/FCS/TakeOff/P Height'
 * '<S234>' : 'mainV03_56/FCS/TakeOff/PID Controller'
 * '<S235>' : 'mainV03_56/FCS/TakeOff/PID Controller1'
 * '<S236>' : 'mainV03_56/FCS/TakeOff/PID Controller2'
 * '<S237>' : 'mainV03_56/FCS/TakeOff/PID Controller3'
 * '<S238>' : 'mainV03_56/FCS/TakeOff/PID Controller4'
 * '<S239>' : 'mainV03_56/FCS/TakeOff/PID Controller5'
 * '<S240>' : 'mainV03_56/FCS/TakeOff/Rotation Angles to Direction Cosine Matrix'
 * '<S241>' : 'mainV03_56/FCS/TakeOff/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S242>' : 'mainV03_56/LIBIS Plant/Aerodynamics'
 * '<S243>' : 'mainV03_56/LIBIS Plant/Custom Variable Mass 6DOF (Euler Angles) (propio)'
 * '<S244>' : 'mainV03_56/LIBIS Plant/Flat Earth to LLA'
 * '<S245>' : 'mainV03_56/LIBIS Plant/Incidence, Sideslip'
 * '<S246>' : 'mainV03_56/LIBIS Plant/Inertia Block'
 * '<S247>' : 'mainV03_56/LIBIS Plant/Propulsion'
 * '<S248>' : 'mainV03_56/LIBIS Plant/Subsystem'
 * '<S249>' : 'mainV03_56/LIBIS Plant/TransferFunction'
 * '<S250>' : 'mainV03_56/LIBIS Plant/TransferFunctions'
 * '<S251>' : 'mainV03_56/LIBIS Plant/alpha_dot, beta_dot'
 * '<S252>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Aerodynamic Forces and Moments '
 * '<S253>' : 'mainV03_56/LIBIS Plant/Aerodynamics/CD'
 * '<S254>' : 'mainV03_56/LIBIS Plant/Aerodynamics/CL'
 * '<S255>' : 'mainV03_56/LIBIS Plant/Aerodynamics/CY'
 * '<S256>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Cl'
 * '<S257>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Cm'
 * '<S258>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Cn'
 * '<S259>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Create Transformation'
 * '<S260>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Dynamic Pressure'
 * '<S261>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Incidence, Sideslip, & Airspeed'
 * '<S262>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Aerodynamic Forces and Moments /3x3 Cross Product'
 * '<S263>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Aerodynamic Forces and Moments /CG-CP Transformation'
 * '<S264>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Aerodynamic Forces and Moments /Force Transformation'
 * '<S265>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Aerodynamic Forces and Moments /Moment Transformation'
 * '<S266>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Aerodynamic Forces and Moments /3x3 Cross Product/Subsystem'
 * '<S267>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Aerodynamic Forces and Moments /3x3 Cross Product/Subsystem1'
 * '<S268>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Aerodynamic Forces and Moments /CG-CP Transformation/Create Transformation'
 * '<S269>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Aerodynamic Forces and Moments /CG-CP Transformation/Incidence, Sideslip, & Airspeed'
 * '<S270>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Aerodynamic Forces and Moments /CG-CP Transformation/Create Transformation/Create 3x3 Matrix'
 * '<S271>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Aerodynamic Forces and Moments /CG-CP Transformation/Incidence, Sideslip, & Airspeed/Subsystem'
 * '<S272>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Aerodynamic Forces and Moments /CG-CP Transformation/Incidence, Sideslip, & Airspeed/Subsystem1'
 * '<S273>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Aerodynamic Forces and Moments /CG-CP Transformation/Incidence, Sideslip, & Airspeed/dot'
 * '<S274>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Aerodynamic Forces and Moments /Force Transformation/Create Transformation'
 * '<S275>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Aerodynamic Forces and Moments /Force Transformation/Incidence, Sideslip, & Airspeed'
 * '<S276>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Aerodynamic Forces and Moments /Force Transformation/Create Transformation/Create 3x3 Matrix'
 * '<S277>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Aerodynamic Forces and Moments /Force Transformation/Incidence, Sideslip, & Airspeed/Subsystem'
 * '<S278>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Aerodynamic Forces and Moments /Force Transformation/Incidence, Sideslip, & Airspeed/Subsystem1'
 * '<S279>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Aerodynamic Forces and Moments /Force Transformation/Incidence, Sideslip, & Airspeed/dot'
 * '<S280>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Aerodynamic Forces and Moments /Moment Transformation/Create Transformation'
 * '<S281>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Aerodynamic Forces and Moments /Moment Transformation/Incidence, Sideslip, & Airspeed'
 * '<S282>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Aerodynamic Forces and Moments /Moment Transformation/Create Transformation/Create 3x3 Matrix'
 * '<S283>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Aerodynamic Forces and Moments /Moment Transformation/Incidence, Sideslip, & Airspeed/Subsystem'
 * '<S284>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Aerodynamic Forces and Moments /Moment Transformation/Incidence, Sideslip, & Airspeed/Subsystem1'
 * '<S285>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Aerodynamic Forces and Moments /Moment Transformation/Incidence, Sideslip, & Airspeed/dot'
 * '<S286>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Create Transformation/Create 3x3 Matrix'
 * '<S287>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Dynamic Pressure/dot'
 * '<S288>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Incidence, Sideslip, & Airspeed/Subsystem'
 * '<S289>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Incidence, Sideslip, & Airspeed/Subsystem1'
 * '<S290>' : 'mainV03_56/LIBIS Plant/Aerodynamics/Incidence, Sideslip, & Airspeed/dot'
 * '<S291>' : 'mainV03_56/LIBIS Plant/Custom Variable Mass 6DOF (Euler Angles) (propio)/Calculate DCM & Euler Angles'
 * '<S292>' : 'mainV03_56/LIBIS Plant/Custom Variable Mass 6DOF (Euler Angles) (propio)/Calculate omega_dot'
 * '<S293>' : 'mainV03_56/LIBIS Plant/Custom Variable Mass 6DOF (Euler Angles) (propio)/Determine Force,  Mass & Inertia'
 * '<S294>' : 'mainV03_56/LIBIS Plant/Custom Variable Mass 6DOF (Euler Angles) (propio)/Vbxw'
 * '<S295>' : 'mainV03_56/LIBIS Plant/Custom Variable Mass 6DOF (Euler Angles) (propio)/Velocity Conversion'
 * '<S296>' : 'mainV03_56/LIBIS Plant/Custom Variable Mass 6DOF (Euler Angles) (propio)/Velocity Conversion1'
 * '<S297>' : 'mainV03_56/LIBIS Plant/Custom Variable Mass 6DOF (Euler Angles) (propio)/Velocity Conversion2'
 * '<S298>' : 'mainV03_56/LIBIS Plant/Custom Variable Mass 6DOF (Euler Angles) (propio)/transform to Inertial axes '
 * '<S299>' : 'mainV03_56/LIBIS Plant/Custom Variable Mass 6DOF (Euler Angles) (propio)/Calculate DCM & Euler Angles/Rotation Angles to Direction Cosine Matrix'
 * '<S300>' : 'mainV03_56/LIBIS Plant/Custom Variable Mass 6DOF (Euler Angles) (propio)/Calculate DCM & Euler Angles/phidot thetadot psidot'
 * '<S301>' : 'mainV03_56/LIBIS Plant/Custom Variable Mass 6DOF (Euler Angles) (propio)/Calculate DCM & Euler Angles/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S302>' : 'mainV03_56/LIBIS Plant/Custom Variable Mass 6DOF (Euler Angles) (propio)/Calculate omega_dot/I x w'
 * '<S303>' : 'mainV03_56/LIBIS Plant/Custom Variable Mass 6DOF (Euler Angles) (propio)/Calculate omega_dot/I x w1'
 * '<S304>' : 'mainV03_56/LIBIS Plant/Custom Variable Mass 6DOF (Euler Angles) (propio)/Calculate omega_dot/wx(Iw)'
 * '<S305>' : 'mainV03_56/LIBIS Plant/Custom Variable Mass 6DOF (Euler Angles) (propio)/Calculate omega_dot/wx(Iw)/Subsystem'
 * '<S306>' : 'mainV03_56/LIBIS Plant/Custom Variable Mass 6DOF (Euler Angles) (propio)/Calculate omega_dot/wx(Iw)/Subsystem1'
 * '<S307>' : 'mainV03_56/LIBIS Plant/Custom Variable Mass 6DOF (Euler Angles) (propio)/Determine Force,  Mass & Inertia/Mass input//output  momentum'
 * '<S308>' : 'mainV03_56/LIBIS Plant/Custom Variable Mass 6DOF (Euler Angles) (propio)/Determine Force,  Mass & Inertia/Mass input//output  momentum/For Each Subsystem'
 * '<S309>' : 'mainV03_56/LIBIS Plant/Custom Variable Mass 6DOF (Euler Angles) (propio)/Vbxw/Subsystem'
 * '<S310>' : 'mainV03_56/LIBIS Plant/Custom Variable Mass 6DOF (Euler Angles) (propio)/Vbxw/Subsystem1'
 * '<S311>' : 'mainV03_56/LIBIS Plant/Flat Earth to LLA/LatLong wrap'
 * '<S312>' : 'mainV03_56/LIBIS Plant/Flat Earth to LLA/LatLong wrap1'
 * '<S313>' : 'mainV03_56/LIBIS Plant/Flat Earth to LLA/LongLat_offset'
 * '<S314>' : 'mainV03_56/LIBIS Plant/Flat Earth to LLA/pos_deg'
 * '<S315>' : 'mainV03_56/LIBIS Plant/Flat Earth to LLA/LatLong wrap/Latitude Wrap 90'
 * '<S316>' : 'mainV03_56/LIBIS Plant/Flat Earth to LLA/LatLong wrap/Wrap Longitude'
 * '<S317>' : 'mainV03_56/LIBIS Plant/Flat Earth to LLA/LatLong wrap/Latitude Wrap 90/Compare To Constant'
 * '<S318>' : 'mainV03_56/LIBIS Plant/Flat Earth to LLA/LatLong wrap/Latitude Wrap 90/Wrap Angle 180'
 * '<S319>' : 'mainV03_56/LIBIS Plant/Flat Earth to LLA/LatLong wrap/Latitude Wrap 90/Wrap Angle 180/Compare To Constant'
 * '<S320>' : 'mainV03_56/LIBIS Plant/Flat Earth to LLA/LatLong wrap/Wrap Longitude/Compare To Constant'
 * '<S321>' : 'mainV03_56/LIBIS Plant/Flat Earth to LLA/LatLong wrap1/Latitude Wrap 90'
 * '<S322>' : 'mainV03_56/LIBIS Plant/Flat Earth to LLA/LatLong wrap1/Wrap Longitude'
 * '<S323>' : 'mainV03_56/LIBIS Plant/Flat Earth to LLA/LatLong wrap1/Latitude Wrap 90/Compare To Constant'
 * '<S324>' : 'mainV03_56/LIBIS Plant/Flat Earth to LLA/LatLong wrap1/Latitude Wrap 90/Wrap Angle 180'
 * '<S325>' : 'mainV03_56/LIBIS Plant/Flat Earth to LLA/LatLong wrap1/Latitude Wrap 90/Wrap Angle 180/Compare To Constant'
 * '<S326>' : 'mainV03_56/LIBIS Plant/Flat Earth to LLA/LatLong wrap1/Wrap Longitude/Compare To Constant'
 * '<S327>' : 'mainV03_56/LIBIS Plant/Flat Earth to LLA/LongLat_offset/Angle Conversion2'
 * '<S328>' : 'mainV03_56/LIBIS Plant/Flat Earth to LLA/LongLat_offset/Find Radian//Distance'
 * '<S329>' : 'mainV03_56/LIBIS Plant/Flat Earth to LLA/LongLat_offset/Find Radian//Distance/Angle Conversion2'
 * '<S330>' : 'mainV03_56/LIBIS Plant/Flat Earth to LLA/LongLat_offset/Find Radian//Distance/denom'
 * '<S331>' : 'mainV03_56/LIBIS Plant/Flat Earth to LLA/LongLat_offset/Find Radian//Distance/e'
 * '<S332>' : 'mainV03_56/LIBIS Plant/Flat Earth to LLA/LongLat_offset/Find Radian//Distance/e^4'
 * '<S333>' : 'mainV03_56/LIBIS Plant/Incidence, Sideslip/Subsystem'
 * '<S334>' : 'mainV03_56/LIBIS Plant/Incidence, Sideslip/Subsystem1'
 * '<S335>' : 'mainV03_56/LIBIS Plant/Incidence, Sideslip/dot'
 * '<S336>' : 'mainV03_56/LIBIS Plant/Inertia Block/Estimate Center of Gravity'
 * '<S337>' : 'mainV03_56/LIBIS Plant/Inertia Block/Symmetric Inertia Tensor'
 * '<S338>' : 'mainV03_56/LIBIS Plant/Inertia Block/Symmetric Inertia Tensor1'
 * '<S339>' : 'mainV03_56/LIBIS Plant/Inertia Block/Estimate Center of Gravity/Interpolate CG'
 * '<S340>' : 'mainV03_56/LIBIS Plant/Inertia Block/Estimate Center of Gravity/Interpolate CG/Matrix interpolation'
 * '<S341>' : 'mainV03_56/LIBIS Plant/Inertia Block/Symmetric Inertia Tensor/Create 3x3 Matrix'
 * '<S342>' : 'mainV03_56/LIBIS Plant/Inertia Block/Symmetric Inertia Tensor1/Create 3x3 Matrix'
 * '<S343>' : 'mainV03_56/LIBIS Plant/Propulsion/8S 25C Battery'
 * '<S344>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_1'
 * '<S345>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_2'
 * '<S346>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_3'
 * '<S347>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_4'
 * '<S348>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_5'
 * '<S349>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_1/mappedMotor_1'
 * '<S350>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_1/noMappedMotor_1'
 * '<S351>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_1/mappedMotor_1/3x3 Cross Product'
 * '<S352>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_1/mappedMotor_1/ESC_1'
 * '<S353>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_1/mappedMotor_1/Motor and Propeller Coupling_1'
 * '<S354>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_1/mappedMotor_1/Torque'
 * '<S355>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_1/mappedMotor_1/3x3 Cross Product/Subsystem'
 * '<S356>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_1/mappedMotor_1/3x3 Cross Product/Subsystem1'
 * '<S357>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_1/noMappedMotor_1/ESC_1'
 * '<S358>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_1/noMappedMotor_1/Motor_1'
 * '<S359>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_1/noMappedMotor_1/Propeller_1'
 * '<S360>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_1/noMappedMotor_1/Torque'
 * '<S361>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_1/noMappedMotor_1/Motor_1/ActiveMotor1'
 * '<S362>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_1/noMappedMotor_1/Motor_1/noActiveMotor1'
 * '<S363>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_2/mappedMotor_2'
 * '<S364>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_2/noMappedMotor_2'
 * '<S365>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_2/mappedMotor_2/3x3 Cross Product'
 * '<S366>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_2/mappedMotor_2/ESC_2'
 * '<S367>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_2/mappedMotor_2/Motor and Propeller Coupling_2'
 * '<S368>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_2/mappedMotor_2/Torque'
 * '<S369>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_2/mappedMotor_2/3x3 Cross Product/Subsystem'
 * '<S370>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_2/mappedMotor_2/3x3 Cross Product/Subsystem1'
 * '<S371>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_2/noMappedMotor_2/ESC_2'
 * '<S372>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_2/noMappedMotor_2/Motor_2'
 * '<S373>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_2/noMappedMotor_2/Propeller_2'
 * '<S374>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_2/noMappedMotor_2/Torque'
 * '<S375>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_2/noMappedMotor_2/Motor_2/ActiveMotor2'
 * '<S376>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_2/noMappedMotor_2/Motor_2/noActiveMotor2'
 * '<S377>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_3/mappedMotor_3'
 * '<S378>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_3/noMappedMotor_3'
 * '<S379>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_3/mappedMotor_3/3x3 Cross Product'
 * '<S380>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_3/mappedMotor_3/ESC_3'
 * '<S381>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_3/mappedMotor_3/Motor and Propeller Coupling_3'
 * '<S382>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_3/mappedMotor_3/Torque'
 * '<S383>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_3/mappedMotor_3/3x3 Cross Product/Subsystem'
 * '<S384>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_3/mappedMotor_3/3x3 Cross Product/Subsystem1'
 * '<S385>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_3/noMappedMotor_3/ESC_3'
 * '<S386>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_3/noMappedMotor_3/Motor_3'
 * '<S387>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_3/noMappedMotor_3/Propeller_3'
 * '<S388>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_3/noMappedMotor_3/Torque'
 * '<S389>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_3/noMappedMotor_3/Motor_3/ActiveMotor3'
 * '<S390>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_3/noMappedMotor_3/Motor_3/noActiveMotor3'
 * '<S391>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_4/mappedMotor_4'
 * '<S392>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_4/noMappedMotor_4'
 * '<S393>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_4/mappedMotor_4/3x3 Cross Product'
 * '<S394>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_4/mappedMotor_4/ESC_4'
 * '<S395>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_4/mappedMotor_4/Motor and Propeller Coupling_4'
 * '<S396>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_4/mappedMotor_4/Torque'
 * '<S397>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_4/mappedMotor_4/3x3 Cross Product/Subsystem'
 * '<S398>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_4/mappedMotor_4/3x3 Cross Product/Subsystem1'
 * '<S399>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_4/noMappedMotor_4/ESC_4'
 * '<S400>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_4/noMappedMotor_4/Motor_4'
 * '<S401>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_4/noMappedMotor_4/Propeller_4'
 * '<S402>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_4/noMappedMotor_4/Torque'
 * '<S403>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_4/noMappedMotor_4/Motor_4/ActiveMotor4'
 * '<S404>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_4/noMappedMotor_4/Motor_4/noActiveMotor4'
 * '<S405>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_5/mappedMotor_5'
 * '<S406>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_5/noMappedMotor_5'
 * '<S407>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_5/mappedMotor_5/3x3 Cross Product'
 * '<S408>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_5/mappedMotor_5/ESC_5'
 * '<S409>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_5/mappedMotor_5/Motor and Propeller Coupling_5'
 * '<S410>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_5/mappedMotor_5/Torque'
 * '<S411>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_5/mappedMotor_5/3x3 Cross Product/Subsystem'
 * '<S412>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_5/mappedMotor_5/3x3 Cross Product/Subsystem1'
 * '<S413>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_5/noMappedMotor_5/ESC_5'
 * '<S414>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_5/noMappedMotor_5/Motor_5'
 * '<S415>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_5/noMappedMotor_5/Propeller_5'
 * '<S416>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_5/noMappedMotor_5/Torque'
 * '<S417>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_5/noMappedMotor_5/Motor_5/ActiveMotor5'
 * '<S418>' : 'mainV03_56/LIBIS Plant/Propulsion/Variant_Motor_5/noMappedMotor_5/Motor_5/noActiveMotor5'
 * '<S419>' : 'mainV03_56/LIBIS Plant/Subsystem/Direction Cosine Matrix Body to Wind'
 * '<S420>' : 'mainV03_56/LIBIS Plant/Subsystem/calc alpha_dot'
 * '<S421>' : 'mainV03_56/LIBIS Plant/Subsystem/Direction Cosine Matrix Body to Wind/A11'
 * '<S422>' : 'mainV03_56/LIBIS Plant/Subsystem/Direction Cosine Matrix Body to Wind/A12'
 * '<S423>' : 'mainV03_56/LIBIS Plant/Subsystem/Direction Cosine Matrix Body to Wind/A13'
 * '<S424>' : 'mainV03_56/LIBIS Plant/Subsystem/Direction Cosine Matrix Body to Wind/A21'
 * '<S425>' : 'mainV03_56/LIBIS Plant/Subsystem/Direction Cosine Matrix Body to Wind/A22'
 * '<S426>' : 'mainV03_56/LIBIS Plant/Subsystem/Direction Cosine Matrix Body to Wind/A23'
 * '<S427>' : 'mainV03_56/LIBIS Plant/Subsystem/Direction Cosine Matrix Body to Wind/A31'
 * '<S428>' : 'mainV03_56/LIBIS Plant/Subsystem/Direction Cosine Matrix Body to Wind/A32'
 * '<S429>' : 'mainV03_56/LIBIS Plant/Subsystem/Direction Cosine Matrix Body to Wind/A33'
 * '<S430>' : 'mainV03_56/LIBIS Plant/Subsystem/calc alpha_dot/Direction Cosine Matrix Body to Wind'
 * '<S431>' : 'mainV03_56/LIBIS Plant/Subsystem/calc alpha_dot/Direction Cosine Matrix Body to Wind/A11'
 * '<S432>' : 'mainV03_56/LIBIS Plant/Subsystem/calc alpha_dot/Direction Cosine Matrix Body to Wind/A12'
 * '<S433>' : 'mainV03_56/LIBIS Plant/Subsystem/calc alpha_dot/Direction Cosine Matrix Body to Wind/A13'
 * '<S434>' : 'mainV03_56/LIBIS Plant/Subsystem/calc alpha_dot/Direction Cosine Matrix Body to Wind/A21'
 * '<S435>' : 'mainV03_56/LIBIS Plant/Subsystem/calc alpha_dot/Direction Cosine Matrix Body to Wind/A22'
 * '<S436>' : 'mainV03_56/LIBIS Plant/Subsystem/calc alpha_dot/Direction Cosine Matrix Body to Wind/A23'
 * '<S437>' : 'mainV03_56/LIBIS Plant/Subsystem/calc alpha_dot/Direction Cosine Matrix Body to Wind/A31'
 * '<S438>' : 'mainV03_56/LIBIS Plant/Subsystem/calc alpha_dot/Direction Cosine Matrix Body to Wind/A32'
 * '<S439>' : 'mainV03_56/LIBIS Plant/Subsystem/calc alpha_dot/Direction Cosine Matrix Body to Wind/A33'
 * '<S440>' : 'mainV03_56/LIBIS Plant/alpha_dot, beta_dot/Subsystem'
 * '<S441>' : 'mainV03_56/LIBIS Plant/alpha_dot, beta_dot/Subsystem/Direction Cosine Matrix Body to Wind'
 * '<S442>' : 'mainV03_56/LIBIS Plant/alpha_dot, beta_dot/Subsystem/calc alpha_dot'
 * '<S443>' : 'mainV03_56/LIBIS Plant/alpha_dot, beta_dot/Subsystem/Direction Cosine Matrix Body to Wind/A11'
 * '<S444>' : 'mainV03_56/LIBIS Plant/alpha_dot, beta_dot/Subsystem/Direction Cosine Matrix Body to Wind/A12'
 * '<S445>' : 'mainV03_56/LIBIS Plant/alpha_dot, beta_dot/Subsystem/Direction Cosine Matrix Body to Wind/A13'
 * '<S446>' : 'mainV03_56/LIBIS Plant/alpha_dot, beta_dot/Subsystem/Direction Cosine Matrix Body to Wind/A21'
 * '<S447>' : 'mainV03_56/LIBIS Plant/alpha_dot, beta_dot/Subsystem/Direction Cosine Matrix Body to Wind/A22'
 * '<S448>' : 'mainV03_56/LIBIS Plant/alpha_dot, beta_dot/Subsystem/Direction Cosine Matrix Body to Wind/A23'
 * '<S449>' : 'mainV03_56/LIBIS Plant/alpha_dot, beta_dot/Subsystem/Direction Cosine Matrix Body to Wind/A31'
 * '<S450>' : 'mainV03_56/LIBIS Plant/alpha_dot, beta_dot/Subsystem/Direction Cosine Matrix Body to Wind/A32'
 * '<S451>' : 'mainV03_56/LIBIS Plant/alpha_dot, beta_dot/Subsystem/Direction Cosine Matrix Body to Wind/A33'
 * '<S452>' : 'mainV03_56/LIBIS Plant/alpha_dot, beta_dot/Subsystem/calc alpha_dot/Direction Cosine Matrix Body to Wind'
 * '<S453>' : 'mainV03_56/LIBIS Plant/alpha_dot, beta_dot/Subsystem/calc alpha_dot/Direction Cosine Matrix Body to Wind/A11'
 * '<S454>' : 'mainV03_56/LIBIS Plant/alpha_dot, beta_dot/Subsystem/calc alpha_dot/Direction Cosine Matrix Body to Wind/A12'
 * '<S455>' : 'mainV03_56/LIBIS Plant/alpha_dot, beta_dot/Subsystem/calc alpha_dot/Direction Cosine Matrix Body to Wind/A13'
 * '<S456>' : 'mainV03_56/LIBIS Plant/alpha_dot, beta_dot/Subsystem/calc alpha_dot/Direction Cosine Matrix Body to Wind/A21'
 * '<S457>' : 'mainV03_56/LIBIS Plant/alpha_dot, beta_dot/Subsystem/calc alpha_dot/Direction Cosine Matrix Body to Wind/A22'
 * '<S458>' : 'mainV03_56/LIBIS Plant/alpha_dot, beta_dot/Subsystem/calc alpha_dot/Direction Cosine Matrix Body to Wind/A23'
 * '<S459>' : 'mainV03_56/LIBIS Plant/alpha_dot, beta_dot/Subsystem/calc alpha_dot/Direction Cosine Matrix Body to Wind/A31'
 * '<S460>' : 'mainV03_56/LIBIS Plant/alpha_dot, beta_dot/Subsystem/calc alpha_dot/Direction Cosine Matrix Body to Wind/A32'
 * '<S461>' : 'mainV03_56/LIBIS Plant/alpha_dot, beta_dot/Subsystem/calc alpha_dot/Direction Cosine Matrix Body to Wind/A33'
 * '<S462>' : 'mainV03_56/Pilot Commands/Buttons to Switches'
 * '<S463>' : 'mainV03_56/Pilot Commands/Signal Builder'
 * '<S464>' : 'mainV03_56/Pilot Commands/Buttons to Switches/Clock'
 * '<S465>' : 'mainV03_56/Pilot Commands/Buttons to Switches/J-K Flip-Flop 1'
 * '<S466>' : 'mainV03_56/Pilot Commands/Buttons to Switches/J-K Flip-Flop 10'
 * '<S467>' : 'mainV03_56/Pilot Commands/Buttons to Switches/J-K Flip-Flop 11'
 * '<S468>' : 'mainV03_56/Pilot Commands/Buttons to Switches/J-K Flip-Flop 12'
 * '<S469>' : 'mainV03_56/Pilot Commands/Buttons to Switches/J-K Flip-Flop 2'
 * '<S470>' : 'mainV03_56/Pilot Commands/Buttons to Switches/J-K Flip-Flop 3'
 * '<S471>' : 'mainV03_56/Pilot Commands/Buttons to Switches/J-K Flip-Flop 4'
 * '<S472>' : 'mainV03_56/Pilot Commands/Buttons to Switches/J-K Flip-Flop 5'
 * '<S473>' : 'mainV03_56/Pilot Commands/Buttons to Switches/J-K Flip-Flop 6'
 * '<S474>' : 'mainV03_56/Pilot Commands/Buttons to Switches/J-K Flip-Flop 7'
 * '<S475>' : 'mainV03_56/Pilot Commands/Buttons to Switches/J-K Flip-Flop 8'
 * '<S476>' : 'mainV03_56/Pilot Commands/Buttons to Switches/J-K Flip-Flop 9'
 * '<S477>' : 'mainV03_56/Pilot Commands/Buttons to Switches/J-K Flip-Flop 1/J-K Flip-Flop'
 * '<S478>' : 'mainV03_56/Pilot Commands/Buttons to Switches/J-K Flip-Flop 10/J-K Flip-Flop'
 * '<S479>' : 'mainV03_56/Pilot Commands/Buttons to Switches/J-K Flip-Flop 11/J-K Flip-Flop'
 * '<S480>' : 'mainV03_56/Pilot Commands/Buttons to Switches/J-K Flip-Flop 12/J-K Flip-Flop'
 * '<S481>' : 'mainV03_56/Pilot Commands/Buttons to Switches/J-K Flip-Flop 2/J-K Flip-Flop'
 * '<S482>' : 'mainV03_56/Pilot Commands/Buttons to Switches/J-K Flip-Flop 3/J-K Flip-Flop'
 * '<S483>' : 'mainV03_56/Pilot Commands/Buttons to Switches/J-K Flip-Flop 4/J-K Flip-Flop'
 * '<S484>' : 'mainV03_56/Pilot Commands/Buttons to Switches/J-K Flip-Flop 5/J-K Flip-Flop'
 * '<S485>' : 'mainV03_56/Pilot Commands/Buttons to Switches/J-K Flip-Flop 6/J-K Flip-Flop'
 * '<S486>' : 'mainV03_56/Pilot Commands/Buttons to Switches/J-K Flip-Flop 7/J-K Flip-Flop'
 * '<S487>' : 'mainV03_56/Pilot Commands/Buttons to Switches/J-K Flip-Flop 8/J-K Flip-Flop'
 * '<S488>' : 'mainV03_56/Pilot Commands/Buttons to Switches/J-K Flip-Flop 9/J-K Flip-Flop'
 * '<S489>' : 'mainV03_56/Sensors/Feedthrough Sensors'
 * '<S490>' : 'mainV03_56/Sensors/Noisy Sensors'
 * '<S491>' : 'mainV03_56/Sensors/Feedthrough Sensors/CPU'
 * '<S492>' : 'mainV03_56/Sensors/Feedthrough Sensors/GPS'
 * '<S493>' : 'mainV03_56/Sensors/Feedthrough Sensors/IMU'
 * '<S494>' : 'mainV03_56/Sensors/Feedthrough Sensors/Pitot'
 * '<S495>' : 'mainV03_56/Sensors/Feedthrough Sensors/Veleta'
 * '<S496>' : 'mainV03_56/Sensors/Feedthrough Sensors/IMU/Calculate DCM & Euler Angles'
 * '<S497>' : 'mainV03_56/Sensors/Noisy Sensors/CPU (Feedthough)'
 * '<S498>' : 'mainV03_56/Sensors/Noisy Sensors/GPS (Feedthough)'
 * '<S499>' : 'mainV03_56/Sensors/Noisy Sensors/IMU'
 * '<S500>' : 'mainV03_56/Sensors/Noisy Sensors/Pitot (Feedthough)'
 * '<S501>' : 'mainV03_56/Sensors/Noisy Sensors/Veleta (Feedthough)'
 * '<S502>' : 'mainV03_56/Sensors/Noisy Sensors/GPS (Feedthough)/LLA to Flat Earth'
 * '<S503>' : 'mainV03_56/Sensors/Noisy Sensors/GPS (Feedthough)/LLA to Flat Earth/LatLong wrap'
 * '<S504>' : 'mainV03_56/Sensors/Noisy Sensors/GPS (Feedthough)/LLA to Flat Earth/LatLong wrap LL0'
 * '<S505>' : 'mainV03_56/Sensors/Noisy Sensors/GPS (Feedthough)/LLA to Flat Earth/Subsystem'
 * '<S506>' : 'mainV03_56/Sensors/Noisy Sensors/GPS (Feedthough)/LLA to Flat Earth/pos_rad'
 * '<S507>' : 'mainV03_56/Sensors/Noisy Sensors/GPS (Feedthough)/LLA to Flat Earth/LatLong wrap/Latitude Wrap 90'
 * '<S508>' : 'mainV03_56/Sensors/Noisy Sensors/GPS (Feedthough)/LLA to Flat Earth/LatLong wrap/Wrap Longitude'
 * '<S509>' : 'mainV03_56/Sensors/Noisy Sensors/GPS (Feedthough)/LLA to Flat Earth/LatLong wrap/Latitude Wrap 90/Compare To Constant'
 * '<S510>' : 'mainV03_56/Sensors/Noisy Sensors/GPS (Feedthough)/LLA to Flat Earth/LatLong wrap/Latitude Wrap 90/Wrap Angle 180'
 * '<S511>' : 'mainV03_56/Sensors/Noisy Sensors/GPS (Feedthough)/LLA to Flat Earth/LatLong wrap/Latitude Wrap 90/Wrap Angle 180/Compare To Constant'
 * '<S512>' : 'mainV03_56/Sensors/Noisy Sensors/GPS (Feedthough)/LLA to Flat Earth/LatLong wrap/Wrap Longitude/Compare To Constant'
 * '<S513>' : 'mainV03_56/Sensors/Noisy Sensors/GPS (Feedthough)/LLA to Flat Earth/LatLong wrap LL0/Latitude Wrap 90'
 * '<S514>' : 'mainV03_56/Sensors/Noisy Sensors/GPS (Feedthough)/LLA to Flat Earth/LatLong wrap LL0/Wrap Longitude'
 * '<S515>' : 'mainV03_56/Sensors/Noisy Sensors/GPS (Feedthough)/LLA to Flat Earth/LatLong wrap LL0/Latitude Wrap 90/Compare To Constant'
 * '<S516>' : 'mainV03_56/Sensors/Noisy Sensors/GPS (Feedthough)/LLA to Flat Earth/LatLong wrap LL0/Latitude Wrap 90/Wrap Angle 180'
 * '<S517>' : 'mainV03_56/Sensors/Noisy Sensors/GPS (Feedthough)/LLA to Flat Earth/LatLong wrap LL0/Latitude Wrap 90/Wrap Angle 180/Compare To Constant'
 * '<S518>' : 'mainV03_56/Sensors/Noisy Sensors/GPS (Feedthough)/LLA to Flat Earth/LatLong wrap LL0/Wrap Longitude/Compare To Constant'
 * '<S519>' : 'mainV03_56/Sensors/Noisy Sensors/GPS (Feedthough)/LLA to Flat Earth/Subsystem/Angle Conversion2'
 * '<S520>' : 'mainV03_56/Sensors/Noisy Sensors/GPS (Feedthough)/LLA to Flat Earth/Subsystem/Find Radian//Distance'
 * '<S521>' : 'mainV03_56/Sensors/Noisy Sensors/GPS (Feedthough)/LLA to Flat Earth/Subsystem/Find Radian//Distance/Angle Conversion2'
 * '<S522>' : 'mainV03_56/Sensors/Noisy Sensors/GPS (Feedthough)/LLA to Flat Earth/Subsystem/Find Radian//Distance/denom'
 * '<S523>' : 'mainV03_56/Sensors/Noisy Sensors/GPS (Feedthough)/LLA to Flat Earth/Subsystem/Find Radian//Distance/e'
 * '<S524>' : 'mainV03_56/Sensors/Noisy Sensors/GPS (Feedthough)/LLA to Flat Earth/Subsystem/Find Radian//Distance/e^4'
 * '<S525>' : 'mainV03_56/Sensors/Noisy Sensors/IMU/Calculate DCM & Euler Angles'
 * '<S526>' : 'mainV03_56/Sensors/Noisy Sensors/IMU/Three-axis Inertial Measurement Unit'
 * '<S527>' : 'mainV03_56/Sensors/Noisy Sensors/IMU/Calculate DCM & Euler Angles/Rotation Angles to Direction Cosine Matrix'
 * '<S528>' : 'mainV03_56/Sensors/Noisy Sensors/IMU/Calculate DCM & Euler Angles/phidot thetadot psidot'
 * '<S529>' : 'mainV03_56/Sensors/Noisy Sensors/IMU/Calculate DCM & Euler Angles/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S530>' : 'mainV03_56/Sensors/Noisy Sensors/IMU/Three-axis Inertial Measurement Unit/Acceleration Conversion'
 * '<S531>' : 'mainV03_56/Sensors/Noisy Sensors/IMU/Three-axis Inertial Measurement Unit/Three-axis Accelerometer'
 * '<S532>' : 'mainV03_56/Sensors/Noisy Sensors/IMU/Three-axis Inertial Measurement Unit/Three-axis Gyroscope'
 * '<S533>' : 'mainV03_56/Sensors/Noisy Sensors/IMU/Three-axis Inertial Measurement Unit/Three-axis Accelerometer/Dynamics'
 * '<S534>' : 'mainV03_56/Sensors/Noisy Sensors/IMU/Three-axis Inertial Measurement Unit/Three-axis Accelerometer/Random bias'
 * '<S535>' : 'mainV03_56/Sensors/Noisy Sensors/IMU/Three-axis Inertial Measurement Unit/Three-axis Accelerometer/w x (w x d)'
 * '<S536>' : 'mainV03_56/Sensors/Noisy Sensors/IMU/Three-axis Inertial Measurement Unit/Three-axis Accelerometer/wdot x d'
 * '<S537>' : 'mainV03_56/Sensors/Noisy Sensors/IMU/Three-axis Inertial Measurement Unit/Three-axis Accelerometer/Dynamics/No Dynamics'
 * '<S538>' : 'mainV03_56/Sensors/Noisy Sensors/IMU/Three-axis Inertial Measurement Unit/Three-axis Accelerometer/Dynamics/Second-order Dynamics'
 * '<S539>' : 'mainV03_56/Sensors/Noisy Sensors/IMU/Three-axis Inertial Measurement Unit/Three-axis Accelerometer/w x (w x d)/w x (w x d)'
 * '<S540>' : 'mainV03_56/Sensors/Noisy Sensors/IMU/Three-axis Inertial Measurement Unit/Three-axis Accelerometer/w x (w x d)/w x d'
 * '<S541>' : 'mainV03_56/Sensors/Noisy Sensors/IMU/Three-axis Inertial Measurement Unit/Three-axis Accelerometer/w x (w x d)/w x (w x d)/Subsystem'
 * '<S542>' : 'mainV03_56/Sensors/Noisy Sensors/IMU/Three-axis Inertial Measurement Unit/Three-axis Accelerometer/w x (w x d)/w x (w x d)/Subsystem1'
 * '<S543>' : 'mainV03_56/Sensors/Noisy Sensors/IMU/Three-axis Inertial Measurement Unit/Three-axis Accelerometer/w x (w x d)/w x d/Subsystem'
 * '<S544>' : 'mainV03_56/Sensors/Noisy Sensors/IMU/Three-axis Inertial Measurement Unit/Three-axis Accelerometer/w x (w x d)/w x d/Subsystem1'
 * '<S545>' : 'mainV03_56/Sensors/Noisy Sensors/IMU/Three-axis Inertial Measurement Unit/Three-axis Accelerometer/wdot x d/Subsystem'
 * '<S546>' : 'mainV03_56/Sensors/Noisy Sensors/IMU/Three-axis Inertial Measurement Unit/Three-axis Accelerometer/wdot x d/Subsystem1'
 * '<S547>' : 'mainV03_56/Sensors/Noisy Sensors/IMU/Three-axis Inertial Measurement Unit/Three-axis Gyroscope/Dynamics'
 * '<S548>' : 'mainV03_56/Sensors/Noisy Sensors/IMU/Three-axis Inertial Measurement Unit/Three-axis Gyroscope/Random bias'
 * '<S549>' : 'mainV03_56/Sensors/Noisy Sensors/IMU/Three-axis Inertial Measurement Unit/Three-axis Gyroscope/Dynamics/No Dynamics'
 * '<S550>' : 'mainV03_56/Sensors/Noisy Sensors/IMU/Three-axis Inertial Measurement Unit/Three-axis Gyroscope/Dynamics/Second-order Dynamics'
 * '<S551>' : 'mainV03_56/Visualization/Axes to VR Axes'
 * '<S552>' : 'mainV03_56/Visualization/Rotation Angles to Direction Cosine Matrix'
 * '<S553>' : 'mainV03_56/Visualization/Rotation Matrix to VRML Rotation'
 * '<S554>' : 'mainV03_56/Visualization/Subsystem'
 * '<S555>' : 'mainV03_56/Visualization/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S556>' : 'mainV03_56/Visualization/Rotation Matrix to VRML Rotation/General case'
 * '<S557>' : 'mainV03_56/Visualization/Rotation Matrix to VRML Rotation/Phi == 0'
 * '<S558>' : 'mainV03_56/Visualization/Rotation Matrix to VRML Rotation/Phi == pi'
 * '<S559>' : 'mainV03_56/Visualization/Rotation Matrix to VRML Rotation/Phi == pi/Logic for flipping axis signs'
 * '<S560>' : 'mainV03_56/Visualization/Rotation Matrix to VRML Rotation/Phi == pi/Logic for flipping axis signs/Compare To Zero'
 * '<S561>' : 'mainV03_56/Visualization/Subsystem/Angle Conversion'
 * '<S562>' : 'mainV03_56/Visualization/Subsystem/Angle Conversion1'
 * '<S563>' : 'mainV03_56/Visualization/Subsystem/Flightpath Angle'
 * '<S564>' : 'mainV03_56/Visualization/Subsystem/Length Conversion'
 * '<S565>' : 'mainV03_56/Visualization/Subsystem/Velocity Conversion'
 * '<S566>' : 'mainV03_56/Visualization/Subsystem/Velocity Conversion1'
 */
#endif                                 /* RTW_HEADER_mainV03_56_h_ */

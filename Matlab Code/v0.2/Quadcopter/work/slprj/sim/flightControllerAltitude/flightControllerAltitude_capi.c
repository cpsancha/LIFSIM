#include "__cf_flightControllerAltitude.h"
#include "rtw_capi.h"
#ifdef HOST_CAPI_BUILD
#include "flightControllerAltitude_capi_host.h"
#define sizeof(s) ((size_t)(0xFFFF))
#undef rt_offsetof
#define rt_offsetof(s,el) ((uint16_T)(0xFFFF))
#define TARGET_CONST
#define TARGET_STRING(s) (s)    
#else
#include "builtin_typeid_types.h"
#include "flightControllerAltitude.h"
#include "flightControllerAltitude_capi.h"
#include "flightControllerAltitude_private.h"
#ifdef LIGHT_WEIGHT_CAPI
#define TARGET_CONST                  
#define TARGET_STRING(s)               (NULL)                    
#else
#define TARGET_CONST                   const
#define TARGET_STRING(s)               (s)
#endif
#endif
static rtwCAPI_Signals rtBlockSignals [ ] = { { 0 , 0 , ( NULL ) , ( NULL ) ,
0 , 0 , 0 , 0 , 0 } } ; static rtwCAPI_States rtBlockStates [ ] = { { 0 , 0 ,
TARGET_STRING ( "flightControllerAltitude/Controller/PID Altitude/Filter" ) ,
TARGET_STRING ( "" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 , 1 , - 1 , 0
} , { 1 , 1 , TARGET_STRING (
"flightControllerAltitude/Controller/PID Altitude/Integrator" ) ,
TARGET_STRING ( "" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 , 1 , - 1 , 0
} , { 2 , 2 , TARGET_STRING (
"flightControllerAltitude/Controller/PID Roll/Filter" ) , TARGET_STRING ( ""
) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 , 1 , - 1 , 0 } , { 3 , 3 ,
TARGET_STRING ( "flightControllerAltitude/Controller/PID Roll/Integrator" ) ,
TARGET_STRING ( "" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 , 1 , - 1 , 0
} , { 4 , 4 , TARGET_STRING (
"flightControllerAltitude/Controller/PID pitch/Filter" ) , TARGET_STRING ( ""
) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 , 1 , - 1 , 0 } , { 5 , 5 ,
TARGET_STRING ( "flightControllerAltitude/Controller/PID pitch/Integrator" )
, TARGET_STRING ( "" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 , 1 , - 1 ,
0 } , { 0 , - 1 , ( NULL ) , ( NULL ) , ( NULL ) , 0 , 0 , 0 , 0 , 0 , 0 , -
1 , 0 } } ;
#ifndef HOST_CAPI_BUILD
static void flightControllerAltitude_InitializeDataAddr ( void * dataAddr [ ]
, h3hamlavmu * localX ) { dataAddr [ 0 ] = ( void * ) ( & localX ->
bts1x2lksb ) ; dataAddr [ 1 ] = ( void * ) ( & localX -> lt4hgqreiu ) ;
dataAddr [ 2 ] = ( void * ) ( & localX -> j3arznukek ) ; dataAddr [ 3 ] = (
void * ) ( & localX -> kykkwvysrp ) ; dataAddr [ 4 ] = ( void * ) ( & localX
-> o2w4enabdh ) ; dataAddr [ 5 ] = ( void * ) ( & localX -> j4u0puqoim ) ; }
#endif
#ifndef HOST_CAPI_BUILD
static void flightControllerAltitude_InitializeVarDimsAddr ( int32_T *
vardimsAddr [ ] ) { vardimsAddr [ 0 ] = ( NULL ) ; }
#endif
static TARGET_CONST rtwCAPI_DataTypeMap rtDataTypeMap [ ] = { { "double" ,
"real_T" , 0 , 0 , sizeof ( real_T ) , SS_DOUBLE , 0 , 0 } } ;
#ifdef HOST_CAPI_BUILD
#undef sizeof
#endif
static TARGET_CONST rtwCAPI_ElementMap rtElementMap [ ] = { { ( NULL ) , 0 ,
0 , 0 , 0 } , } ; static rtwCAPI_DimensionMap rtDimensionMap [ ] = { {
rtwCAPI_SCALAR , 0 , 2 , 0 } } ; static uint_T rtDimensionArray [ ] = { 1 , 1
} ; static const real_T rtcapiStoredFloats [ ] = { 0.0 } ; static
rtwCAPI_FixPtMap rtFixPtMap [ ] = { { ( NULL ) , ( NULL ) ,
rtwCAPI_FIX_RESERVED , 0 , 0 , 0 } , } ; static rtwCAPI_SampleTimeMap
rtSampleTimeMap [ ] = { { ( const void * ) & rtcapiStoredFloats [ 0 ] , (
const void * ) & rtcapiStoredFloats [ 0 ] , 0 , 0 } } ; static int_T
rtContextSystems [ 2 ] ; static rtwCAPI_LoggingMetaInfo loggingMetaInfo [ ] =
{ { 0 , 0 , "" , 0 } } ; static rtwCAPI_ModelMapLoggingStaticInfo
mmiStaticInfoLogging = { 2 , rtContextSystems , loggingMetaInfo , 0 , NULL ,
{ 0 , NULL , NULL } , 0 , ( NULL ) } ; static rtwCAPI_ModelMappingStaticInfo
mmiStatic = { { rtBlockSignals , 0 , ( NULL ) , 0 , ( NULL ) , 0 } , { ( NULL
) , 0 , ( NULL ) , 0 } , { rtBlockStates , 6 } , { rtDataTypeMap ,
rtDimensionMap , rtFixPtMap , rtElementMap , rtSampleTimeMap ,
rtDimensionArray } , "float" , { 1221716897U , 1938829669U , 3831392184U ,
1199966055U } , & mmiStaticInfoLogging , 0 , 0 } ; const
rtwCAPI_ModelMappingStaticInfo * flightControllerAltitude_GetCAPIStaticMap (
) { return & mmiStatic ; }
#ifndef HOST_CAPI_BUILD
static void flightControllerAltitude_InitializeSystemRan ( f0y0h2ugli * const
mpsvw5q1fx , sysRanDType * systemRan [ ] , int_T systemTid [ ] , void *
rootSysRanPtr , int rootTid ) { UNUSED_PARAMETER ( mpsvw5q1fx ) ; systemRan [
0 ] = ( sysRanDType * ) rootSysRanPtr ; systemRan [ 1 ] = ( NULL ) ;
systemTid [ 1 ] = hm1aavqkmk [ 0 ] ; systemTid [ 0 ] = rootTid ;
rtContextSystems [ 0 ] = 0 ; rtContextSystems [ 1 ] = 0 ; }
#endif
#ifndef HOST_CAPI_BUILD
void flightControllerAltitude_InitializeDataMapInfo ( f0y0h2ugli * const
mpsvw5q1fx , h3hamlavmu * localX , void * sysRanPtr , int contextTid ) {
rtwCAPI_SetVersion ( mpsvw5q1fx -> DataMapInfo . mmi , 1 ) ;
rtwCAPI_SetStaticMap ( mpsvw5q1fx -> DataMapInfo . mmi , & mmiStatic ) ;
rtwCAPI_SetLoggingStaticMap ( mpsvw5q1fx -> DataMapInfo . mmi , &
mmiStaticInfoLogging ) ; flightControllerAltitude_InitializeDataAddr (
mpsvw5q1fx -> DataMapInfo . dataAddress , localX ) ;
rtwCAPI_SetDataAddressMap ( mpsvw5q1fx -> DataMapInfo . mmi , mpsvw5q1fx ->
DataMapInfo . dataAddress ) ; flightControllerAltitude_InitializeVarDimsAddr
( mpsvw5q1fx -> DataMapInfo . vardimsAddress ) ; rtwCAPI_SetVarDimsAddressMap
( mpsvw5q1fx -> DataMapInfo . mmi , mpsvw5q1fx -> DataMapInfo .
vardimsAddress ) ; rtwCAPI_SetPath ( mpsvw5q1fx -> DataMapInfo . mmi , ( NULL
) ) ; rtwCAPI_SetFullPath ( mpsvw5q1fx -> DataMapInfo . mmi , ( NULL ) ) ;
rtwCAPI_SetInstanceLoggingInfo ( mpsvw5q1fx -> DataMapInfo . mmi , &
mpsvw5q1fx -> DataMapInfo . mmiLogInstanceInfo ) ; rtwCAPI_SetChildMMIArray (
mpsvw5q1fx -> DataMapInfo . mmi , ( NULL ) ) ; rtwCAPI_SetChildMMIArrayLen (
mpsvw5q1fx -> DataMapInfo . mmi , 0 ) ;
flightControllerAltitude_InitializeSystemRan ( mpsvw5q1fx , mpsvw5q1fx ->
DataMapInfo . systemRan , mpsvw5q1fx -> DataMapInfo . systemTid , sysRanPtr ,
contextTid ) ; rtwCAPI_SetSystemRan ( mpsvw5q1fx -> DataMapInfo . mmi ,
mpsvw5q1fx -> DataMapInfo . systemRan ) ; rtwCAPI_SetSystemTid ( mpsvw5q1fx
-> DataMapInfo . mmi , mpsvw5q1fx -> DataMapInfo . systemTid ) ;
rtwCAPI_SetGlobalTIDMap ( mpsvw5q1fx -> DataMapInfo . mmi , & hm1aavqkmk [ 0
] ) ; }
#else
#ifdef __cplusplus
extern "C" {
#endif
void flightControllerAltitude_host_InitializeDataMapInfo (
flightControllerAltitude_host_DataMapInfo_T * dataMap , const char * path ) {
rtwCAPI_SetVersion ( dataMap -> mmi , 1 ) ; rtwCAPI_SetStaticMap ( dataMap ->
mmi , & mmiStatic ) ; rtwCAPI_SetDataAddressMap ( dataMap -> mmi , NULL ) ;
rtwCAPI_SetVarDimsAddressMap ( dataMap -> mmi , NULL ) ; rtwCAPI_SetPath (
dataMap -> mmi , path ) ; rtwCAPI_SetFullPath ( dataMap -> mmi , NULL ) ;
rtwCAPI_SetChildMMIArray ( dataMap -> mmi , ( NULL ) ) ;
rtwCAPI_SetChildMMIArrayLen ( dataMap -> mmi , 0 ) ; }
#ifdef __cplusplus
}
#endif
#endif

#include "__cf_flightControllerAltitude.h"
#ifndef RTW_HEADER_flightControllerAltitude_h_
#define RTW_HEADER_flightControllerAltitude_h_
#include <string.h>
#include <stddef.h>
#include "rtw_modelmap.h"
#ifndef flightControllerAltitude_COMMON_INCLUDES_
#define flightControllerAltitude_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "simstruc.h"
#include "fixedpoint.h"
#endif
#include "flightControllerAltitude_types.h"
#include "multiword_types.h"
#include "rt_nonfinite.h"
typedef struct { real_T eanqaqsd1j ; real_T mrbvuxzexc ; real_T hcmj4t5rih ;
real_T ilpj4pbyou ; real_T iusm0z5qjj ; real_T oulygiyoye ; } nuontokzu0 ;
typedef struct { real_T bts1x2lksb ; real_T lt4hgqreiu ; real_T j3arznukek ;
real_T kykkwvysrp ; real_T o2w4enabdh ; real_T j4u0puqoim ; } h3hamlavmu ;
typedef struct { real_T bts1x2lksb ; real_T lt4hgqreiu ; real_T j3arznukek ;
real_T kykkwvysrp ; real_T o2w4enabdh ; real_T j4u0puqoim ; } lrg2bkmtsu ;
typedef struct { boolean_T bts1x2lksb ; boolean_T lt4hgqreiu ; boolean_T
j3arznukek ; boolean_T kykkwvysrp ; boolean_T o2w4enabdh ; boolean_T
j4u0puqoim ; } oircpoxb2p ; struct hayblwnioy_ { real_T P_1 ; real_T P_2 ;
real_T P_3 ; real_T P_4 ; real_T P_5 ; real_T P_6 ; real_T P_7 ; real_T P_8 ;
real_T P_9 ; real_T P_10 ; real_T P_11 ; real_T P_12 ; real_T P_13 ; real_T
P_14 ; real_T P_15 ; real_T P_16 ; real_T P_17 ; real_T P_18 ; real_T P_19 ;
real_T P_20 ; real_T P_21 ; real_T P_22 ; real_T P_23 ; real_T P_24 ; real_T
P_25 ; real_T P_26 ; real_T P_27 ; } ; struct bwi5gh53uv { struct
SimStruct_tag * _mdlRefSfcnS ; struct { rtwCAPI_ModelMappingInfo mmi ;
rtwCAPI_ModelMapLoggingInstanceInfo mmiLogInstanceInfo ; void * dataAddress [
6 ] ; int32_T * vardimsAddress [ 6 ] ; sysRanDType * systemRan [ 2 ] ; int_T
systemTid [ 2 ] ; } DataMapInfo ; } ; typedef struct { nuontokzu0 rtb ;
f0y0h2ugli rtm ; } h5j0o4j2xxd ; extern real_T rtP_altitudeHoldValue ; extern
void d3ewoot30z ( SimStruct * _mdlRefSfcnS , int_T mdlref_TID0 , int_T
mdlref_TID1 , int_T mdlref_TID2 , f0y0h2ugli * const mpsvw5q1fx , nuontokzu0
* localB , h3hamlavmu * localX , void * sysRanPtr , int contextTid ,
rtwCAPI_ModelMappingInfo * rt_ParentMMI , const char_T * rt_ChildPath , int_T
rt_ChildMMIIdx , int_T rt_CSTATEIdx ) ; extern void
mr_flightControllerAltitude_MdlInfoRegFcn ( SimStruct * mdlRefSfcnS , char_T
* modelName , int_T * retVal ) ; extern mxArray *
mr_flightControllerAltitude_GetDWork ( const h5j0o4j2xxd * mdlrefDW ) ;
extern void mr_flightControllerAltitude_SetDWork ( h5j0o4j2xxd * mdlrefDW ,
const mxArray * ssDW ) ; extern void
mr_flightControllerAltitude_RegisterSimStateChecksum ( SimStruct * S ) ;
extern mxArray * mr_flightControllerAltitude_GetSimStateDisallowedBlocks ( )
; extern const rtwCAPI_ModelMappingStaticInfo *
flightControllerAltitude_GetCAPIStaticMap ( void ) ; extern void egikcehx0x (
h3hamlavmu * localX ) ; extern void nehev0e1d1 ( h3hamlavmu * localX ) ;
extern void okpmulkuep ( nuontokzu0 * localB , lrg2bkmtsu * localXdot ) ;
extern void ogpba5xu4y ( void ) ; extern void ogpba5xu4yTID2 ( void ) ;
extern void flightControllerAltitude ( const real_T * e0pjqdljq0 , const
real_T * fkzuw0qg0a , const real_T * faeskpdsrx , const SensorsBus *
hp0zhdhr0h , real_T c5ra12ogi1 [ 4 ] , nuontokzu0 * localB , h3hamlavmu *
localX ) ; extern void flightControllerAltitudeTID2 ( void ) ;
#endif

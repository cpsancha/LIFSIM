#include "__cf_flightControllerAltitude.h"
#ifndef RTW_HEADER_flightControllerAltitude_cap_host_h_
#define RTW_HEADER_flightControllerAltitude_cap_host_h_
#ifdef HOST_CAPI_BUILD
#include "rtw_capi.h"
#include "rtw_modelmap.h"
typedef struct { rtwCAPI_ModelMappingInfo mmi ; }
flightControllerAltitude_host_DataMapInfo_T ;
#ifdef __cplusplus
extern "C" {
#endif
void flightControllerAltitude_host_InitializeDataMapInfo (
flightControllerAltitude_host_DataMapInfo_T * dataMap , const char * path ) ;
#ifdef __cplusplus
}
#endif
#endif
#endif

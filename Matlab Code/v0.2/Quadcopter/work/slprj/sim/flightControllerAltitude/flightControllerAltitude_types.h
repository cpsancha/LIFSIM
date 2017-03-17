#include "__cf_flightControllerAltitude.h"
#ifndef RTW_HEADER_flightControllerAltitude_types_h_
#define RTW_HEADER_flightControllerAltitude_types_h_
#include "rtwtypes.h"
#include "builtin_typeid_types.h"
#include "multiword_types.h"
#ifndef DEFINED_TYPEDEF_FOR_CommandBus_
#define DEFINED_TYPEDEF_FOR_CommandBus_
typedef struct { real_T roll_cmd ; real_T pitch_cmd ; real_T yaw_cmd ; real_T
throttle_cmd ; } CommandBus ;
#endif
#ifndef DEFINED_TYPEDEF_FOR_SensorsBus_
#define DEFINED_TYPEDEF_FOR_SensorsBus_
typedef struct { real_T AccelMeas_body [ 3 ] ; real_T OmegaMeas_body [ 3 ] ;
real_T LLA [ 3 ] ; } SensorsBus ;
#endif
#ifndef DEFINED_TYPEDEF_FOR_struct_VTjambXRg1S46wFaJTouEB_
#define DEFINED_TYPEDEF_FOR_struct_VTjambXRg1S46wFaJTouEB_
typedef struct { real_T dPitch ; } struct_VTjambXRg1S46wFaJTouEB ;
#endif
#ifndef DEFINED_TYPEDEF_FOR_struct_mjd6v4gaUxSn9lufjJDe7_
#define DEFINED_TYPEDEF_FOR_struct_mjd6v4gaUxSn9lufjJDe7_
typedef struct { real_T dRoll ; } struct_mjd6v4gaUxSn9lufjJDe7 ;
#endif
#ifndef DEFINED_TYPEDEF_FOR_struct_NWzN4KJIFUNsigjCGU5IzE_
#define DEFINED_TYPEDEF_FOR_struct_NWzN4KJIFUNsigjCGU5IzE_
typedef struct { real_T dYaw ; } struct_NWzN4KJIFUNsigjCGU5IzE ;
#endif
#ifndef DEFINED_TYPEDEF_FOR_struct_5mlylREK2pOgdEHL9dH1yE_
#define DEFINED_TYPEDEF_FOR_struct_5mlylREK2pOgdEHL9dH1yE_
typedef struct { real_T dThrottle ; real_T bias ; }
struct_5mlylREK2pOgdEHL9dH1yE ;
#endif
#ifndef DEFINED_TYPEDEF_FOR_struct_h6h7nIv6T9OCWPyBw2OrpG_
#define DEFINED_TYPEDEF_FOR_struct_h6h7nIv6T9OCWPyBw2OrpG_
typedef struct { struct_VTjambXRg1S46wFaJTouEB Pitch ;
struct_mjd6v4gaUxSn9lufjJDe7 Roll ; struct_NWzN4KJIFUNsigjCGU5IzE Yaw ;
struct_5mlylREK2pOgdEHL9dH1yE Throttle ; } struct_h6h7nIv6T9OCWPyBw2OrpG ;
#endif
#ifndef DEFINED_TYPEDEF_FOR_struct_H2zzByA25P0LxlJ9Q30baC_
#define DEFINED_TYPEDEF_FOR_struct_H2zzByA25P0LxlJ9Q30baC_
typedef struct { struct_h6h7nIv6T9OCWPyBw2OrpG QuadMix ; }
struct_H2zzByA25P0LxlJ9Q30baC ;
#endif
typedef struct hayblwnioy_ hayblwnioy ; typedef struct bwi5gh53uv f0y0h2ugli
;
#endif

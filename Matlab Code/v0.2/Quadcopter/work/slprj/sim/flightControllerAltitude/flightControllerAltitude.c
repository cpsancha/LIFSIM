#include "__cf_flightControllerAltitude.h"
#include "flightControllerAltitude_capi.h"
#include "flightControllerAltitude.h"
#include "flightControllerAltitude_private.h"
int_T hm1aavqkmk [ 3 ] ; static RegMdlInfo rtMdlInfo_flightControllerAltitude
[ 51 ] = { { "h5j0o4j2xxd" , MDL_INFO_NAME_MDLREF_DWORK , 0 , - 1 , ( void *
) "flightControllerAltitude" } , { "lc5221ncc4" ,
MDL_INFO_ID_GLOBAL_RTW_CONSTRUCT , 0 , - 1 , ( void * )
"flightControllerAltitude" } , { "fr0mrnpg0a" ,
MDL_INFO_ID_GLOBAL_RTW_CONSTRUCT , 0 , - 1 , ( void * )
"flightControllerAltitude" } , { "oircpoxb2p" ,
MDL_INFO_ID_GLOBAL_RTW_CONSTRUCT , 0 , - 1 , ( void * )
"flightControllerAltitude" } , { "lrg2bkmtsu" ,
MDL_INFO_ID_GLOBAL_RTW_CONSTRUCT , 0 , - 1 , ( void * )
"flightControllerAltitude" } , { "h3hamlavmu" ,
MDL_INFO_ID_GLOBAL_RTW_CONSTRUCT , 0 , - 1 , ( void * )
"flightControllerAltitude" } , { "gsbb1pl3kh" ,
MDL_INFO_ID_GLOBAL_RTW_CONSTRUCT , 0 , - 1 , ( void * )
"flightControllerAltitude" } , { "fji4tw4wcn" ,
MDL_INFO_ID_GLOBAL_RTW_CONSTRUCT , 0 , - 1 , ( void * )
"flightControllerAltitude" } , { "gjlnvd1nla" ,
MDL_INFO_ID_GLOBAL_RTW_CONSTRUCT , 0 , - 1 , ( void * )
"flightControllerAltitude" } , { "flc4cy2ebj" ,
MDL_INFO_ID_GLOBAL_RTW_CONSTRUCT , 0 , - 1 , ( void * )
"flightControllerAltitude" } , { "jxleqitvce" ,
MDL_INFO_ID_GLOBAL_RTW_CONSTRUCT , 0 , - 1 , ( void * )
"flightControllerAltitude" } , { "nuontokzu0" ,
MDL_INFO_ID_GLOBAL_RTW_CONSTRUCT , 0 , - 1 , ( void * )
"flightControllerAltitude" } , { "kgbwhvjl02" ,
MDL_INFO_ID_GLOBAL_RTW_CONSTRUCT , 0 , - 1 , ( void * )
"flightControllerAltitude" } , { "d3dw32bw2r" ,
MDL_INFO_ID_GLOBAL_RTW_CONSTRUCT , 0 , - 1 , ( void * )
"flightControllerAltitude" } , { "jvyqxiv531" ,
MDL_INFO_ID_GLOBAL_RTW_CONSTRUCT , 0 , - 1 , ( void * )
"flightControllerAltitude" } , { "okpmulkuep" ,
MDL_INFO_ID_GLOBAL_RTW_CONSTRUCT , 0 , - 1 , ( void * )
"flightControllerAltitude" } , { "ogpba5xu4y" ,
MDL_INFO_ID_GLOBAL_RTW_CONSTRUCT , 0 , - 1 , ( void * )
"flightControllerAltitude" } , { "nehev0e1d1" ,
MDL_INFO_ID_GLOBAL_RTW_CONSTRUCT , 0 , - 1 , ( void * )
"flightControllerAltitude" } , { "egikcehx0x" ,
MDL_INFO_ID_GLOBAL_RTW_CONSTRUCT , 0 , - 1 , ( void * )
"flightControllerAltitude" } , { "d3ewoot30z" ,
MDL_INFO_ID_GLOBAL_RTW_CONSTRUCT , 0 , - 1 , ( void * )
"flightControllerAltitude" } , { "m00d4lafok" ,
MDL_INFO_ID_GLOBAL_RTW_CONSTRUCT , 0 , - 1 , ( void * )
"flightControllerAltitude" } , { "l20o2fijqg" ,
MDL_INFO_ID_GLOBAL_RTW_CONSTRUCT , 0 , - 1 , ( void * )
"flightControllerAltitude" } , { "flightControllerAltitude" ,
MDL_INFO_ID_GLOBAL_RTW_CONSTRUCT , 0 , 0 , ( NULL ) } , { "bwi5gh53uv" ,
MDL_INFO_ID_GLOBAL_RTW_CONSTRUCT , 0 , - 1 , ( void * )
"flightControllerAltitude" } , { "f0y0h2ugli" ,
MDL_INFO_ID_GLOBAL_RTW_CONSTRUCT , 0 , - 1 , ( void * )
"flightControllerAltitude" } , { "oa0ck3e3n5x" ,
MDL_INFO_ID_GLOBAL_RTW_CONSTRUCT , 0 , - 1 , ( void * )
"flightControllerAltitude" } , { "hm1aavqkmk" ,
MDL_INFO_ID_GLOBAL_RTW_CONSTRUCT , 0 , - 1 , ( void * )
"flightControllerAltitude" } , { "drr5aslbpzr" ,
MDL_INFO_ID_GLOBAL_RTW_CONSTRUCT , 0 , - 1 , ( void * )
"flightControllerAltitude" } , { "oa0ck3e3n5" ,
MDL_INFO_ID_GLOBAL_RTW_CONSTRUCT , 0 , - 1 , ( void * )
"flightControllerAltitude" } , { "struct_H2zzByA25P0LxlJ9Q30baC" ,
MDL_INFO_ID_DATA_TYPE , 0 , - 1 , ( NULL ) } , {
"struct_h6h7nIv6T9OCWPyBw2OrpG" , MDL_INFO_ID_DATA_TYPE , 0 , - 1 , ( NULL )
} , { "struct_5mlylREK2pOgdEHL9dH1yE" , MDL_INFO_ID_DATA_TYPE , 0 , - 1 , (
NULL ) } , { "struct_NWzN4KJIFUNsigjCGU5IzE" , MDL_INFO_ID_DATA_TYPE , 0 , -
1 , ( NULL ) } , { "struct_mjd6v4gaUxSn9lufjJDe7" , MDL_INFO_ID_DATA_TYPE , 0
, - 1 , ( NULL ) } , { "struct_VTjambXRg1S46wFaJTouEB" ,
MDL_INFO_ID_DATA_TYPE , 0 , - 1 , ( NULL ) } , { "SensorsBus" ,
MDL_INFO_ID_DATA_TYPE , 0 , - 1 , ( NULL ) } , { "CommandBus" ,
MDL_INFO_ID_DATA_TYPE , 0 , - 1 , ( NULL ) } , {
"mr_flightControllerAltitude_GetSimStateDisallowedBlocks" ,
MDL_INFO_ID_MODEL_FCN_NAME , 0 , - 1 , ( void * ) "flightControllerAltitude"
} , { "mr_flightControllerAltitude_extractBitFieldFromCellArrayWithOffset" ,
MDL_INFO_ID_MODEL_FCN_NAME , 0 , - 1 , ( void * ) "flightControllerAltitude"
} , { "mr_flightControllerAltitude_cacheBitFieldToCellArrayWithOffset" ,
MDL_INFO_ID_MODEL_FCN_NAME , 0 , - 1 , ( void * ) "flightControllerAltitude"
} , { "mr_flightControllerAltitude_restoreDataFromMxArrayWithOffset" ,
MDL_INFO_ID_MODEL_FCN_NAME , 0 , - 1 , ( void * ) "flightControllerAltitude"
} , { "mr_flightControllerAltitude_cacheDataToMxArrayWithOffset" ,
MDL_INFO_ID_MODEL_FCN_NAME , 0 , - 1 , ( void * ) "flightControllerAltitude"
} , { "mr_flightControllerAltitude_extractBitFieldFromMxArray" ,
MDL_INFO_ID_MODEL_FCN_NAME , 0 , - 1 , ( void * ) "flightControllerAltitude"
} , { "mr_flightControllerAltitude_cacheBitFieldToMxArray" ,
MDL_INFO_ID_MODEL_FCN_NAME , 0 , - 1 , ( void * ) "flightControllerAltitude"
} , { "mr_flightControllerAltitude_restoreDataFromMxArray" ,
MDL_INFO_ID_MODEL_FCN_NAME , 0 , - 1 , ( void * ) "flightControllerAltitude"
} , { "mr_flightControllerAltitude_cacheDataAsMxArray" ,
MDL_INFO_ID_MODEL_FCN_NAME , 0 , - 1 , ( void * ) "flightControllerAltitude"
} , { "mr_flightControllerAltitude_RegisterSimStateChecksum" ,
MDL_INFO_ID_MODEL_FCN_NAME , 0 , - 1 , ( void * ) "flightControllerAltitude"
} , { "mr_flightControllerAltitude_SetDWork" , MDL_INFO_ID_MODEL_FCN_NAME , 0
, - 1 , ( void * ) "flightControllerAltitude" } , {
"mr_flightControllerAltitude_GetDWork" , MDL_INFO_ID_MODEL_FCN_NAME , 0 , - 1
, ( void * ) "flightControllerAltitude" } , { "flightControllerAltitude.h" ,
MDL_INFO_MODEL_FILENAME , 0 , - 1 , ( NULL ) } , {
"flightControllerAltitude.c" , MDL_INFO_MODEL_FILENAME , 0 , - 1 , ( void * )
"flightControllerAltitude" } } ; hayblwnioy hayblwnioyb = { -
2891.38451636519 , 0.00024297273150912 , - 52.8378413871279 , -
0.999560878115567 , 0.000621245528289962 , - 146.118873160313 ,
171.684827700763 , 155.924967958506 , 8.8075624027162 , - 114.050466183378 ,
0.000785869263148391 , - 449.652681173872 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0
, 0.001 , - 1.0 , 0.001 , - 1.0 , 0.005 , 1173.5836838 , 0.001 , - 1.0 , -
1.0 } ; void egikcehx0x ( h3hamlavmu * localX ) { localX -> bts1x2lksb =
hayblwnioyb . P_13 ; localX -> lt4hgqreiu = hayblwnioyb . P_14 ; localX ->
j3arznukek = hayblwnioyb . P_15 ; localX -> kykkwvysrp = hayblwnioyb . P_16 ;
localX -> o2w4enabdh = hayblwnioyb . P_17 ; localX -> j4u0puqoim =
hayblwnioyb . P_18 ; } void nehev0e1d1 ( h3hamlavmu * localX ) { localX ->
bts1x2lksb = hayblwnioyb . P_13 ; localX -> lt4hgqreiu = hayblwnioyb . P_14 ;
localX -> j3arznukek = hayblwnioyb . P_15 ; localX -> kykkwvysrp =
hayblwnioyb . P_16 ; localX -> o2w4enabdh = hayblwnioyb . P_17 ; localX ->
j4u0puqoim = hayblwnioyb . P_18 ; } void flightControllerAltitude ( const
real_T * e0pjqdljq0 , const real_T * fkzuw0qg0a , const real_T * faeskpdsrx ,
const SensorsBus * hp0zhdhr0h , real_T c5ra12ogi1 [ 4 ] , nuontokzu0 * localB
, h3hamlavmu * localX ) { real_T fifd0wo00p ; real_T kerst5n4me ; real_T
nrsy51uchk ; real_T pof2hw011b ; real_T il0khwkukx ; real_T emyup51f5y ;
real_T m20dsgpnsd ; fifd0wo00p = hp0zhdhr0h -> LLA [ 2 ] -
rtP_altitudeHoldValue ; localB -> eanqaqsd1j = ( hayblwnioyb . P_1 *
fifd0wo00p - localX -> bts1x2lksb ) * hayblwnioyb . P_7 ; localB ->
mrbvuxzexc = hayblwnioyb . P_4 * fifd0wo00p ; kerst5n4me = * e0pjqdljq0 -
hp0zhdhr0h -> OmegaMeas_body [ 1 ] ; localB -> hcmj4t5rih = ( hayblwnioyb .
P_2 * kerst5n4me - localX -> j3arznukek ) * hayblwnioyb . P_8 ; localB ->
ilpj4pbyou = hayblwnioyb . P_5 * kerst5n4me ; nrsy51uchk = * fkzuw0qg0a -
hp0zhdhr0h -> OmegaMeas_body [ 0 ] ; localB -> iusm0z5qjj = ( hayblwnioyb .
P_3 * nrsy51uchk - localX -> o2w4enabdh ) * hayblwnioyb . P_9 ; localB ->
oulygiyoye = hayblwnioyb . P_6 * nrsy51uchk ; pof2hw011b = ( ( hayblwnioyb .
P_12 * nrsy51uchk + localX -> j4u0puqoim ) + localB -> iusm0z5qjj ) *
hayblwnioyb . P_19 ; il0khwkukx = ( ( hayblwnioyb . P_11 * kerst5n4me +
localX -> kykkwvysrp ) + localB -> hcmj4t5rih ) * hayblwnioyb . P_21 ;
m20dsgpnsd = ( ( hayblwnioyb . P_10 * fifd0wo00p + localX -> lt4hgqreiu ) +
localB -> eanqaqsd1j ) * hayblwnioyb . P_23 + hayblwnioyb . P_24 ; emyup51f5y
= hayblwnioyb . P_25 * * faeskpdsrx ; c5ra12ogi1 [ 0 ] = ( pof2hw011b +
m20dsgpnsd ) + emyup51f5y ; c5ra12ogi1 [ 1 ] = ( il0khwkukx + m20dsgpnsd ) +
hayblwnioyb . P_26 * emyup51f5y ; c5ra12ogi1 [ 2 ] = ( hayblwnioyb . P_20 *
pof2hw011b + m20dsgpnsd ) + emyup51f5y ; c5ra12ogi1 [ 3 ] = ( hayblwnioyb .
P_22 * il0khwkukx + m20dsgpnsd ) + hayblwnioyb . P_27 * emyup51f5y ; } void
flightControllerAltitudeTID2 ( void ) { } void ogpba5xu4y ( void ) { } void
ogpba5xu4yTID2 ( void ) { } void okpmulkuep ( nuontokzu0 * localB ,
lrg2bkmtsu * localXdot ) { localXdot -> bts1x2lksb = localB -> eanqaqsd1j ;
localXdot -> lt4hgqreiu = localB -> mrbvuxzexc ; localXdot -> j3arznukek =
localB -> hcmj4t5rih ; localXdot -> kykkwvysrp = localB -> ilpj4pbyou ;
localXdot -> o2w4enabdh = localB -> iusm0z5qjj ; localXdot -> j4u0puqoim =
localB -> oulygiyoye ; } void d3ewoot30z ( SimStruct * _mdlRefSfcnS , int_T
mdlref_TID0 , int_T mdlref_TID1 , int_T mdlref_TID2 , f0y0h2ugli * const
mpsvw5q1fx , nuontokzu0 * localB , h3hamlavmu * localX , void * sysRanPtr ,
int contextTid , rtwCAPI_ModelMappingInfo * rt_ParentMMI , const char_T *
rt_ChildPath , int_T rt_ChildMMIIdx , int_T rt_CSTATEIdx ) { rt_InitInfAndNaN
( sizeof ( real_T ) ) ; ( void ) memset ( ( void * ) mpsvw5q1fx , 0 , sizeof
( f0y0h2ugli ) ) ; hm1aavqkmk [ 0 ] = mdlref_TID0 ; hm1aavqkmk [ 1 ] =
mdlref_TID1 ; hm1aavqkmk [ 2 ] = mdlref_TID2 ; mpsvw5q1fx -> _mdlRefSfcnS = (
_mdlRefSfcnS ) ; { localB -> eanqaqsd1j = 0.0 ; localB -> mrbvuxzexc = 0.0 ;
localB -> hcmj4t5rih = 0.0 ; localB -> ilpj4pbyou = 0.0 ; localB ->
iusm0z5qjj = 0.0 ; localB -> oulygiyoye = 0.0 ; }
flightControllerAltitude_InitializeDataMapInfo ( mpsvw5q1fx , localX ,
sysRanPtr , contextTid ) ; if ( ( rt_ParentMMI != ( NULL ) ) && (
rt_ChildPath != ( NULL ) ) ) { rtwCAPI_SetChildMMI ( * rt_ParentMMI ,
rt_ChildMMIIdx , & ( mpsvw5q1fx -> DataMapInfo . mmi ) ) ; rtwCAPI_SetPath (
mpsvw5q1fx -> DataMapInfo . mmi , rt_ChildPath ) ;
rtwCAPI_MMISetContStateStartIndex ( mpsvw5q1fx -> DataMapInfo . mmi ,
rt_CSTATEIdx ) ; } } void mr_flightControllerAltitude_MdlInfoRegFcn (
SimStruct * mdlRefSfcnS , char_T * modelName , int_T * retVal ) { * retVal =
0 ; * retVal = 0 ; ssRegModelRefMdlInfo ( mdlRefSfcnS , modelName ,
rtMdlInfo_flightControllerAltitude , 51 ) ; * retVal = 1 ; } static void
mr_flightControllerAltitude_cacheDataAsMxArray ( mxArray * destArray ,
mwIndex i , int j , const void * srcData , size_t numBytes ) ; static void
mr_flightControllerAltitude_cacheDataAsMxArray ( mxArray * destArray ,
mwIndex i , int j , const void * srcData , size_t numBytes ) { mxArray *
newArray = mxCreateUninitNumericMatrix ( ( size_t ) 1 , numBytes ,
mxUINT8_CLASS , mxREAL ) ; memcpy ( ( uint8_T * ) mxGetData ( newArray ) , (
const uint8_T * ) srcData , numBytes ) ; mxSetFieldByNumber ( destArray , i ,
j , newArray ) ; } static void
mr_flightControllerAltitude_restoreDataFromMxArray ( void * destData , const
mxArray * srcArray , mwIndex i , int j , size_t numBytes ) ; static void
mr_flightControllerAltitude_restoreDataFromMxArray ( void * destData , const
mxArray * srcArray , mwIndex i , int j , size_t numBytes ) { memcpy ( (
uint8_T * ) destData , ( const uint8_T * ) mxGetData ( mxGetFieldByNumber (
srcArray , i , j ) ) , numBytes ) ; } static void
mr_flightControllerAltitude_cacheBitFieldToMxArray ( mxArray * destArray ,
mwIndex i , int j , uint_T bitVal ) ; static void
mr_flightControllerAltitude_cacheBitFieldToMxArray ( mxArray * destArray ,
mwIndex i , int j , uint_T bitVal ) { mxSetFieldByNumber ( destArray , i , j
, mxCreateDoubleScalar ( ( double ) bitVal ) ) ; } static uint_T
mr_flightControllerAltitude_extractBitFieldFromMxArray ( const mxArray *
srcArray , mwIndex i , int j , uint_T numBits ) ; static uint_T
mr_flightControllerAltitude_extractBitFieldFromMxArray ( const mxArray *
srcArray , mwIndex i , int j , uint_T numBits ) { const uint_T varVal = (
uint_T ) mxGetScalar ( mxGetFieldByNumber ( srcArray , i , j ) ) ; return
varVal & ( ( 1u << numBits ) - 1u ) ; } static void
mr_flightControllerAltitude_cacheDataToMxArrayWithOffset ( mxArray *
destArray , mwIndex i , int j , mwIndex offset , const void * srcData ,
size_t numBytes ) ; static void
mr_flightControllerAltitude_cacheDataToMxArrayWithOffset ( mxArray *
destArray , mwIndex i , int j , mwIndex offset , const void * srcData ,
size_t numBytes ) { uint8_T * varData = ( uint8_T * ) mxGetData (
mxGetFieldByNumber ( destArray , i , j ) ) ; memcpy ( ( uint8_T * ) & varData
[ offset * numBytes ] , ( const uint8_T * ) srcData , numBytes ) ; } static
void mr_flightControllerAltitude_restoreDataFromMxArrayWithOffset ( void *
destData , const mxArray * srcArray , mwIndex i , int j , mwIndex offset ,
size_t numBytes ) ; static void
mr_flightControllerAltitude_restoreDataFromMxArrayWithOffset ( void *
destData , const mxArray * srcArray , mwIndex i , int j , mwIndex offset ,
size_t numBytes ) { const uint8_T * varData = ( const uint8_T * ) mxGetData (
mxGetFieldByNumber ( srcArray , i , j ) ) ; memcpy ( ( uint8_T * ) destData ,
( const uint8_T * ) & varData [ offset * numBytes ] , numBytes ) ; } static
void mr_flightControllerAltitude_cacheBitFieldToCellArrayWithOffset ( mxArray
* destArray , mwIndex i , int j , mwIndex offset , uint_T fieldVal ) ; static
void mr_flightControllerAltitude_cacheBitFieldToCellArrayWithOffset ( mxArray
* destArray , mwIndex i , int j , mwIndex offset , uint_T fieldVal ) {
mxSetCell ( mxGetFieldByNumber ( destArray , i , j ) , offset ,
mxCreateDoubleScalar ( ( double ) fieldVal ) ) ; } static uint_T
mr_flightControllerAltitude_extractBitFieldFromCellArrayWithOffset ( const
mxArray * srcArray , mwIndex i , int j , mwIndex offset , uint_T numBits ) ;
static uint_T
mr_flightControllerAltitude_extractBitFieldFromCellArrayWithOffset ( const
mxArray * srcArray , mwIndex i , int j , mwIndex offset , uint_T numBits ) {
const uint_T fieldVal = ( uint_T ) mxGetScalar ( mxGetCell (
mxGetFieldByNumber ( srcArray , i , j ) , offset ) ) ; return fieldVal & ( (
1u << numBits ) - 1u ) ; } mxArray * mr_flightControllerAltitude_GetDWork (
const h5j0o4j2xxd * mdlrefDW ) { static const char * dwFieldNames [ 3 ] = {
"rtb" , "NULL->rtdw" , "NULL->rtzce" , } ; mxArray * ssDW =
mxCreateStructMatrix ( 1 , 1 , 3 , dwFieldNames ) ;
mr_flightControllerAltitude_cacheDataAsMxArray ( ssDW , 0 , 0 , & ( mdlrefDW
-> rtb ) , sizeof ( mdlrefDW -> rtb ) ) ; return ssDW ; } void
mr_flightControllerAltitude_SetDWork ( h5j0o4j2xxd * mdlrefDW , const mxArray
* ssDW ) { mr_flightControllerAltitude_restoreDataFromMxArray ( & ( mdlrefDW
-> rtb ) , ssDW , 0 , 0 , sizeof ( mdlrefDW -> rtb ) ) ; } void
mr_flightControllerAltitude_RegisterSimStateChecksum ( SimStruct * S ) {
const uint32_T chksum [ 4 ] = { 2358562349U , 1494157764U , 86907860U ,
3610805970U , } ; slmrModelRefRegisterSimStateChecksum ( S ,
"flightControllerAltitude" , & chksum [ 0 ] ) ; } mxArray *
mr_flightControllerAltitude_GetSimStateDisallowedBlocks ( ) { return NULL ; }

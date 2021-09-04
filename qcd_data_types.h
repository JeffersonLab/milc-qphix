#ifndef _QCD_DATA_TYPES_H_
#define _QCD_DATA_TYPES_H_

#if QPHIX_PrecisionInt == 1
#ifndef ENABLE_LOW_QPHIX_PrecisionInt
//typedef float fptype;
#undef fptype
#define fptype float
#else
typedef unsigned short half;
#undef fptype
#define fptype half
#endif
#elif QPHIX_PrecisionInt == 2
//typedef double fptype;
#undef fptype
#define fptype double
#else
#error "QPHIX_PrecisionInt Not defined"
#endif

#ifdef COMPRESSED_12
#define COMPRESSED_GAUGES true
#define COMPRESSED_SUFFIX _12
#else
#define COMPRESSED_GAUGES false
#define COMPRESSED_SUFFIX _18
#endif

#ifdef QPHIX_GF
#ifndef GF_PARITY
#define GF_PARITY QPHIX_EVEN
#endif
#endif

template <typename FT, int VL, bool CMPRSD>
class DataTypes {
public:
	typedef FT Spinor[3][4][2][VL];
	typedef FT KS[3][2][VL];

	typedef FT Gauge[8][(CMPRSD? 2 : 3)][3][2][VL];
	typedef FT Gauge18[8][3][3][2][VL];
	typedef FT Gauge12[8][2][3][2][VL];

	typedef FT Hermit[8][8][VL];
	typedef FT HermitHelper[8][7][8][VL];
	typedef FT HermitHelperYZT[6][7][8][VL];

	typedef struct {
		FT diag1[6][VL];
		FT off_diag1[15][2][VL];
		FT diag2[6][VL];
		FT off_diag2[15][2][VL];
	} Clover;

};

#ifndef VECLEN
#error "VECLEN must be defined"
#endif

typedef DataTypes<fptype,VECLEN,COMPRESSED_GAUGES>::Spinor Spinor;
typedef DataTypes<fptype,VECLEN,COMPRESSED_GAUGES>::KS KS;
typedef DataTypes<fptype,VECLEN,COMPRESSED_GAUGES>::Gauge Gauge;
typedef DataTypes<fptype,VECLEN,COMPRESSED_GAUGES>::Gauge18 Gauge18;
typedef DataTypes<fptype,VECLEN,COMPRESSED_GAUGES>::Gauge12 Gauge12;
typedef DataTypes<fptype,VECLEN,COMPRESSED_GAUGES>::Hermit Hermit;
typedef DataTypes<fptype,VECLEN,COMPRESSED_GAUGES>::HermitHelper HermitHelper;
typedef DataTypes<fptype,VECLEN,COMPRESSED_GAUGES>::HermitHelperYZT HermitHelperYZT;
typedef DataTypes<fptype,VECLEN,COMPRESSED_GAUGES>::Clover Clover;

#endif // _QCD_DATA_TYPES_H_

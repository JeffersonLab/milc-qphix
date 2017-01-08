#undef Vxh
#undef Vx
#undef Vy
#undef Vz
#undef Vt
#undef Pxy
#undef Pxyz
#undef BoundTable
#undef NeighTable
#undef BLENGTH
#undef PadBound
#undef PadNeigh
#if QPHIX_PrecisionInt == 1
#define Vxh Vxhf
#define Vx Vxf
#define Vy Vyf
#define Vz Vzf
#define Vt Vtf
#define Pxy Pxyf
#define Pxyz Pxyzf
#define BoundTable BoundTableF
#define NeighTable NeighTableF
#define BLENGTH BLENGTHF
#define PadBound PadBoundF
#define PadNeigh PadNeighF
#define gBar gBarF
#define phaser phaserF
#else
#if QPHIX_PrecisionInt == 2
#define Vxh Vxhd
#define Vx Vxd
#define Vy Vyd
#define Vz Vzd
#define Vt Vtd
#define Pxy Pxyd
#define Pxyz Pxyzd
#define BoundTable BoundTableD
#define NeighTable NeighTableD
#define BLENGTH BLENGTHD
#define PadBound PadBoundD
#define PadNeigh PadNeighD
#define gBar gBarD
#define phaser phaserD
#endif
#endif

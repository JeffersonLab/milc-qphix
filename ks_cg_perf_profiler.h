#ifndef _KS_CG_PERF_PROFILER_H_
#define _KS_CG_PERF_PROFILER_H_

#include "hrtimer.h"   /* High resolution timer */

#ifdef __cplusplus

/* To capture the CG perf metrics */
struct CGPerfData
{
    /* Number of iterations within CG loop. */
    int iters;
    /* The number of elements in the lattice */
    std::size_t latt_size;
    /* Modelled data volumes for the kernels called in the CG routine. */
    std::size_t dslash_bytes;
    std::size_t axpy_bytes;
    std::size_t axpymul_bytes;
    std::size_t axpymul2_bytes;
    std::size_t calc_rsq_bytes;
    std::size_t calc_rsq2_bytes;
    std::size_t calc_rsqmul_bytes;
    std::size_t calc_pkp_bytes;
    std::size_t relative_residue_bytes;
    /* Timing each of the micro-kernels within CG. Time in nanosecs (10^-9) */    
    hrtimer_t dslash_time;
    hrtimer_t axpy_time;
    hrtimer_t calc_rsq_time;
    hrtimer_t calc_rsq2_time;    
//    hrtimer_t calc_rsqmul_time;
    hrtimer_t calc_pkp_time;
    hrtimer_t relative_residue_time;
    
    CGPerfData();
    
    double get_dslash_bw();
    double get_calc_rsq_bw();
    double get_calc_rsq2_bw();
    double get_calc_rsqmul_bw(int);
    double get_calc_pkp_bw();
    double get_relative_residue_bw();
    double get_axpy_bw();
    double get_axpymul_bw(int);
};

inline CGPerfData::CGPerfData()
{
    iters = 0;
    latt_size = Nxh*Ny*Nz*Nt;
    /* Data requirement per iteration of the dslash call*/
    
    /******************* Pre-compute dslash data footprint ********************/
    dslash_bytes = 0;
    /* long links */
#ifndef COMPRESSED_12
    dslash_bytes += 2*latt_size*8*3*3*2*sizeof(fptype);
#else
    dslash_bytes += 2*latt_size*8*2*3*2*sizeof(fptype);
#endif
    /* fat links */
    dslash_bytes += 2*latt_size*8*3*3*2*sizeof(fptype);
    /* ks spinors (assuming we bring in only one spinor per iteration, and reuse 
     * the other seven. */
    dslash_bytes += 2*latt_size*sizeof(KS)/VECLEN; 
    /* output */
    dslash_bytes += 4*latt_size*sizeof(KS)/VECLEN;  

    /***************** pre-compute axpy data footprint ************************/
    axpy_bytes = 3*latt_size*sizeof(KS)/VECLEN; /* y = ax+y */
    
    axpymul_bytes = 4*latt_size*sizeof(KS)/VECLEN; 
    axpymul2_bytes = latt_size*sizeof(KS)/VECLEN;

    /***************** pre-compute calc_rsq data footprint ********************/
    calc_rsq_bytes = 6*latt_size*sizeof(KS)/VECLEN;
    
    /***************** pre-compute calc_rsq2 data footprint *******************/    
    calc_rsq2_bytes = 6*latt_size*sizeof(KS)/VECLEN;
   
    calc_rsqmul_bytes = 3*latt_size*sizeof(KS)/VECLEN;
 
    /***************** pre-compute calc_pkp data footprint ********************/    
    calc_pkp_bytes  = 3*latt_size*sizeof(KS)/VECLEN;
    
    /************* pre-compute relative_residue data footprint ****************/
    relative_residue_bytes = 2*latt_size*sizeof(KS)/VECLEN;
    
    dslash_time    = 0;
    axpy_time      = 0;
    calc_rsq_time  = 0;
    calc_rsq2_time = 0;    
    //calc_rsqmul_time = 0;
    calc_pkp_time  = 0;
    relative_residue_time = 0;
};

double
inline CGPerfData::get_dslash_bw()
{
    return (double)iters*((double)dslash_bytes/(double)dslash_time);
}

double
inline CGPerfData::get_calc_rsq_bw()
{
    return (double)iters*((double)calc_rsq_bytes/(double)calc_rsq_time);
}

double 
inline CGPerfData::get_calc_rsq2_bw()
{
    return (double)iters*((double)calc_rsq2_bytes/(double)calc_rsq2_time);
}

double 
inline CGPerfData::get_calc_rsqmul_bw(int num_offsets)
{
    //return (double)iters*(((num_offsets-1)*(double)calc_rsqmul_bytes+(double)calc_rsq2_bytes)/(double)calc_rsq2_time);
    if(num_offsets>1)
      return (double)iters*(double)calc_rsqmul_bytes/(double)calc_rsq2_time;
    else
      return get_calc_rsq2_bw();
}

double 
inline CGPerfData::get_calc_pkp_bw()
{
    return (double)iters*((double)calc_pkp_bytes/(double)calc_pkp_time);
}

double
inline CGPerfData::get_relative_residue_bw()
{
    return (double)iters*((double)relative_residue_bytes
            /(double)relative_residue_time);
}

double
inline CGPerfData::get_axpy_bw()
{
    return (double)iters*((double)axpy_bytes
            /(double)axpy_time);
}

double
inline CGPerfData::get_axpymul_bw(int num_offsets)
{
    return (double)iters*(((double)axpymul2_bytes+num_offsets*(double)axpymul_bytes)/(double)axpy_time);
}

static void
print_cg_timings(CGPerfData &cgPerfData, int num_offsets=1)
{
		extern int myRank;
    if(myRank == 0) std::cout << "=============================================="
              << "MBENCH CG performance numbers : \n"
              << std::setw(40) << std::left
              << "Total dslash time (s)" << " : "
              << cgPerfData.dslash_time/(1.0e9)
              << std::setw(40)
              << "\nDslash BW (GB/s)" << " : "
              << cgPerfData.get_dslash_bw() 
              << "\n"
              << std::setw(40)
              << "Total calc_rsq time (s)" << " : "
              << cgPerfData.calc_rsq_time/(1.0e9) << "\n"
              << std::setw(40)
              << "Calc RSQ BW (GB/s)" <<  " : "
              << cgPerfData.get_calc_rsq_bw() 
              << "\n"
              << std::setw(40)
              << "Total relative_residue_time (s)" << " : " 
              << cgPerfData.relative_residue_time/(1.0e9) << "\n"
              << std::setw(40)
              << "Relative_residue BW (GB/s)" << " : " << std::setw(10)
              << cgPerfData.get_relative_residue_bw() 
              << "\n"
              << std::setw(40)
              << "Total axpy time (s)" << " : " << std::setw(10)
              << cgPerfData.axpy_time/1.0e9
              << std::setw(40)
              << "\nAxpy BY (GB/s)" << " : " << std::setw(10)
              << cgPerfData.get_axpymul_bw(num_offsets)
              << "\n"
              << std::setw(40)
              << "Total calc_rsq2 time (s)" <<  " : " << std::setw(10)
              << cgPerfData.calc_rsq2_time/(1.0e9)
              << std::setw(40)
              << "\nCalc RSQ2 BW (GB/s)" << " : " << std::setw(10)
              << cgPerfData.get_calc_rsqmul_bw(num_offsets) 
              << std::setw(40)
              << "\nTotal calc_pkp time (s)" << " : " << std::setw(10)
              << cgPerfData.calc_pkp_time/(1.0e9)
              << std::setw(40)
              << "\nCalc PKP BW (GB/s)" << " : " << std::setw(10)
              << cgPerfData.get_calc_pkp_bw() 
              << "\n";
}


#endif

#endif

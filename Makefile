# Makefile to build qphix library

# Must supply these variables:
# ARCH (scalar|knl|hsw|skx)
# PK_CC (choice of C compiler or MPI wrapper)
# PK_CXX (choice of C++ compiler or MPI wrapper)
ARCH ?= scalar

# This file
MAKEFILE=Makefile

# Source directory
BASE_DIR := .

# Output directory
LIB_DIR=${BASE_DIR}/lib

# Installation directory
PREFIX := ../install

# Location of codegen
CODEGEN_PATH := ../milc-qphix-codegen

# Set architecture (knl/avx/sse/scalar)
ARCH:=$(strip $(ARCH))

# Choose to enable debug messages & GDB support (0/1)
ENABLE_DEBUG=0

# Choose to use vtune instrumentation (0/1)
ENABLE_VTUNE=0

# Choose to use MPI (0/1)
ENABLE_MPI=1

TARGET=libqphix

CUSTOM_MAKE=customMake.$(ARCH)
# ifeq ($(ARCH),avx)
#   ifeq ($(AVX2),1)
#     CUSTOM_MAKE=customMake.hsw
#   endif
#   ifeq ($(AVX512),1)
#     # For SKL new custommake has to be created
#     CUSTOM_MAKE=customMake.skl
#   endif
# endif
include $(CUSTOM_MAKE)

ifeq ($(ENABLE_DEBUG),1)
  OPTFLAGS = -g -O0 -DFF_DEBUG -DFF_PROFILE 
else
  OPTFLAGS = -O3
endif

ifeq ($(ENABLE_VTUNE),1)
  OPTFLAGS += -g -DFF_VTUNE
endif

ifeq ($(ENABLE_MPI),1)
  ifneq ($(ASSUME_MULTINODE),1)
    CXX = ${PK_CXX}
  else
    CXX = icpc
  endif
else
  CXX = icpc
endif

CXXFLAGS = -Wall -Wextra $(OPTFLAGS) -qopenmp -fargument-noalias-global -restrict $(PREFFLAGS) -I$(SEPHOME)/include -I./ -I$(CODEGEN_PATH)
CXXFLAGS += -parallel-source-info=2 -debug inline-debug-info -qopt-report=5

LIB_PREFIX := libqphixmilc

# For Xeon SKX
ifeq ($(ARCH),skx)
  VECLENF=16
  VECLEND=8
  
  ifeq ($(AVX512),1) # For SKX
    SUFFIX := _skx
    MPSSFLAGS = -xCORE-AVX512
    yesnolist += AVX512
  else # For KNC
    SUFFIX := _mic
    MPSSFLAGS = -mmic -mGLOB_default_function_attrs="use_gather_scatter_hint=off"
  endif
  LIB_PREFIX := $(LIB_PREFIX)$(SUFFIX)
endif

# For Xeon Phi
ifeq ($(ARCH),knl)
  VECLENF=16
  VECLEND=8
  
  ifeq ($(AVX512),1) # For KNL
    SUFFIX := _avx512
    MPSSFLAGS = -xMIC-AVX512
    yesnolist += AVX512
  else # For KNC
    SUFFIX := _mic
    MPSSFLAGS = -mmic -mGLOB_default_function_attrs="use_gather_scatter_hint=off"
  endif
  LIB_PREFIX := $(LIB_PREFIX)$(SUFFIX)
endif

# For Xeon
ifeq ($(ARCH),avx)
  VECLENF=8
  VECLEND=4
  
  ifeq ($(AVX512),1) # For SKX
    # Suffix avx512 is already taken by KNL
    SUFFIX := _skx
    MPSSFLAGS = -march=core-avx512 -xCORE-AVX512
  endif
  ifeq ($(AVX2),1) # For HSW & BDW
    SUFFIX := _avx2
    MPSSFLAGS = -march=core-avx2 -xCORE-AVX2
  else
    SUFFIX := _avx
    MPSSFLAGS = -mavx
  endif
  LIB_PREFIX := $(LIB_PREFIX)$(SUFFIX)
  MINCT=2
  yesnolist += AVX2
endif

# For older arch
ifeq ($(ARCH),sse)
  VECLENF=4
  VECLEND=2
  
  SUFFIX := _sse
  LIB_PREFIX := $(LIB_PREFIX)$(SUFFIX)
  
  MINCT=2
  MPSSFLAGS = -msse3
endif

# For scalar build
ifeq ($(ARCH),scalar)
  VECLENF=1
  VECLEND=1
  
  SUFFIX := _scalar
  LIB_PREFIX := $(LIB_PREFIX)$(SUFFIX)
endif

ifeq ($(ENABLE_MPI),0)
  ifeq ($(ARCH),knl)
    LIB_PREFIX := $(LIB_PREFIX)_single
    SUFFIX := $(SUFFIX)_single
  endif
  ifeq ($(ARCH),avx)
    LIB_PREFIX := $(LIB_PREFIX)_single
    SUFFIX := $(SUFFIX)_single
  endif
endif
ifeq ($(ENABLE_MPI),1)
  ifeq ($(ARCH),sse)
    LIB_PREFIX := $(LIB_PREFIX)_mpi
    SUFFIX := $(SUFFIX)_mpi
  endif
  ifeq ($(ARCH),scalar)
    LIB_PREFIX := $(LIB_PREFIX)_mpi
    SUFFIX := $(SUFFIX)_mpi
  endif
endif
ifeq ($(ENABLE_DEBUG),1)
    LIB_PREFIX := $(LIB_PREFIX)_dbg
endif
ifeq ($(ENABLE_VTUNE),1)
  LIB_PREFIX := $(LIB_PREFIX)_vt
endif
LIB_FILE := ${LIB_PREFIX}.a
QPHIXLIB := ${LIB_DIR}/${LIB_FILE}


export

yesnolist += ENABLE_STREAMING_STORES
yesnolist += ENABLE_LOW_PRECISION
yesnolist += ENABLE_BARRIER
yesnolist += COMPRESSED_12
yesnolist += USE_WAITANY
yesnolist += ENABLE_MPI
yesnolist += TIME_CG

deflist += MINCT
deflist += BY
deflist += BZ
deflist += BARRIER_FREQ_T
deflist += QPHIX_Precision
deflist += XY_PAD
deflist += XYZ_PAD
deflist += CODEGEN_PATH

DEFS += $(strip $(foreach var, $(yesnolist), $(if $(filter 1, $($(var))), -D$(var))))
DEFS += $(strip $(foreach var, $(deflist), $(if $($(var)), -D$(var)=$($(var)))))

#LIBQPHIXMILC_SRC= ks_long_dslash.cpp ks_d_congrad_fn.cpp  ks_d_congrad_fn_two_src.cpp ks_blas_utils_c.cpp hrtimer.cpp                 
LIBQPHIXMILC_SRC= ks_long_dslash.cpp ks_d_congrad_fn.cpp ks_multicg_offset.cpp ks_blas_utils_c.cpp hrtimer.cpp qphix_f3.cpp qphix_d3.cpp qphix_misc.cpp ks_boundary.cpp
LIBQPHIXMILC_SRC+= gauge_force_imp.cpp gf_boundary.cpp gauge_force_imp_complete_specialization.cpp
LIBQPHIXMILC_SRC+= ff_boundary.cpp fermion_force_complete_specialization.cpp
LIBQPHIXMILC_SRC+= layout.cpp

HEADERS = Barrier.h dslash_mic_complete_specialization.h gauge_force_imp_complete_specialization.h gauge_force_imp.h gf_boundary.h ff_boundary.h fermion_force_complete_specialization.h hrtimer.h ks_blas_utils_c.h ks_boundary.h ks_cg_perf_profiler.h  ks_config.h ks_globals.h phaser.h qcd_data_types.h qphix_d3.h qphix_f3.h qphix.h qphix_internal.h qphix_int.h qphix_types.h stopwatch.h ubench_helper.h utils.h qphix_d3_generic.h qphix_d3_internal.h qphix_f3_generic.h qphix_f3_internal.h ks_long_dslash.h
HEADERS+= fermion_force.h layout.h

all: $(TARGET)

ks_long_dslash.fo: ks_long_dslash.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=1 -DVECLEN=${VECLENF} -c ks_long_dslash.cpp -o ks_long_dslash.fo $(CMLLDFLAGS)

ks_long_dslash.do: ks_long_dslash.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=2 -DVECLEN=${VECLEND} -c ks_long_dslash.cpp -o ks_long_dslash.do $(CMLLDFLAGS)

ks_boundary.fo: ks_boundary.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=1 -DVECLEN=${VECLENF} -c ks_boundary.cpp -o ks_boundary.fo $(CMLLDFLAGS)

ks_boundary.do: ks_boundary.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=2 -DVECLEN=${VECLEND} -c ks_boundary.cpp -o ks_boundary.do $(CMLLDFLAGS)

ks_blas_utils_c.fo: ks_blas_utils_c.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=1 -DVECLEN=${VECLENF} -c ks_blas_utils_c.cpp -o ks_blas_utils_c.fo $(CMLLDFLAGS)

ks_blas_utils_c.do: ks_blas_utils_c.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=2 -DVECLEN=${VECLEND} -c ks_blas_utils_c.cpp -o ks_blas_utils_c.do $(CMLLDFLAGS)

qphix_misc.fo: qphix_misc.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=1 -DVECLEN=${VECLENF} -c qphix_misc.cpp -o qphix_misc.fo $(CMLLDFLAGS)

qphix_misc.do: qphix_misc.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=2 -DVECLEN=${VECLEND} -c qphix_misc.cpp -o qphix_misc.do $(CMLLDFLAGS)

qphix_misc.o: qphix_misc.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=1 -DVECLEN=${VECLENF} -c qphix_misc.cpp -o qphix_misc.o $(CMLLDFLAGS)

ks_d_congrad_fn.fo: ks_d_congrad_fn.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=1 -DVECLEN=${VECLENF} -c ks_d_congrad_fn.cpp -o ks_d_congrad_fn.fo $(CMLLDFLAGS)

ks_d_congrad_fn.do: ks_d_congrad_fn.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=2 -DVECLEN=${VECLEND} -c ks_d_congrad_fn.cpp -o ks_d_congrad_fn.do $(CMLLDFLAGS)

ks_multicg_offset.fo: ks_multicg_offset.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=1 -DVECLEN=${VECLENF} -c ks_multicg_offset.cpp -o ks_multicg_offset.fo $(CMLLDFLAGS)

ks_multicg_offset.do: ks_multicg_offset.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=2 -DVECLEN=${VECLEND} -c ks_multicg_offset.cpp -o ks_multicg_offset.do $(CMLLDFLAGS)

ks_dslash_complete_specialization.fo: ks_dslash_complete_specialization.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=1 -DVECLEN=${VECLENF} -c ks_dslash_complete_specialization.cpp -o ks_dslash_complete_specialization.fo $(CMLLDFLAGS)

ks_dslash_complete_specialization.do: ks_dslash_complete_specialization.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=2 -DVECLEN=${VECLEND} -c ks_dslash_complete_specialization.cpp -o ks_dslash_complete_specialization.do $(CMLLDFLAGS)

#fermion_force_complete_specialization.fo: fermion_force_complete_specialization.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
#	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) -Wno-uninitialized $(DEFS) -DQPHIX_PrecisionInt=1 -DVECLEN=${VECLENF} -c fermion_force_complete_specialization.cpp -o fermion_force_complete_specialization.fo $(CMLLDFLAGS)

fermion_force_complete_specialization.do: fermion_force_complete_specialization.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) -Wno-uninitialized $(DEFS) -DQPHIX_PrecisionInt=2 -DVECLEN=${VECLEND} -c fermion_force_complete_specialization.cpp -o fermion_force_complete_specialization.do $(CMLLDFLAGS)

gauge_force_imp.fo: gauge_force_imp.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=1 -DVECLEN=${VECLENF} -c gauge_force_imp.cpp -o gauge_force_imp.fo $(CMLLDFLAGS)

gauge_force_imp.do: gauge_force_imp.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=2 -DVECLEN=${VECLEND} -c gauge_force_imp.cpp -o gauge_force_imp.do $(CMLLDFLAGS)

gf_boundary.fo: gf_boundary.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=1 -DVECLEN=${VECLENF} -c gf_boundary.cpp -o gf_boundary.fo $(CMLLDFLAGS)

gf_boundary.do: gf_boundary.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=2 -DVECLEN=${VECLEND} -c gf_boundary.cpp -o gf_boundary.do $(CMLLDFLAGS)

#ff_boundary.fo: ff_boundary.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
#	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=1 -DVECLEN=${VECLENF} -c ff_boundary.cpp -o ff_boundary.fo $(CMLLDFLAGS)

ff_boundary.do: ff_boundary.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=2 -DVECLEN=${VECLEND} -c ff_boundary.cpp -o ff_boundary.do $(CMLLDFLAGS)

gauge_force_imp_complete_specialization.fo: gauge_force_imp_complete_specialization.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=1 -DVECLEN=${VECLENF} -c gauge_force_imp_complete_specialization.cpp -o gauge_force_imp_complete_specialization.fo $(CMLLDFLAGS)

gauge_force_imp_complete_specialization.do: gauge_force_imp_complete_specialization.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=2 -DVECLEN=${VECLEND} -c gauge_force_imp_complete_specialization.cpp -o gauge_force_imp_complete_specialization.do $(CMLLDFLAGS)

#layout.fo: layout.c ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
#	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=1 -DVECLEN=${VECLENF} -c layout.c -o layout.fo $(CMLLDFLAGS)

layout.do: layout.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=2 -DVECLEN=${VECLEND} -c layout.cpp -o layout.do $(CMLLDFLAGS) 

#fermion_force.fo: fermion_force.c ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
#	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=1 -DVECLEN=${VECLENF} -c fermion_force.c -o fermion_force.fo $(CMLLDFLAGS)

fermion_force.do: fermion_force.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=2 -DVECLEN=${VECLEND} -DMILC_SM -c fermion_force.cpp -o fermion_force.do $(CMLLDFLAGS) 

#qphix_su3_algebra.fo: qphix_su3_algebra.c ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
#	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=1 -DVECLEN=${VECLENF} -c qphix_su3_algebra.c -o qphix_su3_algebra.fo $(CMLLDFLAGS)

qphix_su3_algebra.do: qphix_su3_algebra.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=2 -DVECLEN=${VECLEND} -DMILC_SM -c qphix_su3_algebra.cpp -o qphix_su3_algebra.do $(CMLLDFLAGS) 

qphix_f3.o: qphix_f3.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=1 -DVECLEN=${VECLENF} -c qphix_f3.cpp -o qphix_f3.o $(CMLLDFLAGS)

qphix_d3.o: qphix_d3.cpp ${HEADERS} $(MAKEFILE) $(CUSTOM_MAKE)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=2 -DVECLEN=${VECLEND} -c qphix_d3.cpp -o qphix_d3.o $(CMLLDFLAGS)

OBJECTS= ks_long_dslash.fo ks_long_dslash.do ks_blas_utils_c.fo ks_blas_utils_c.do qphix_misc.fo qphix_misc.do ks_d_congrad_fn.fo ks_d_congrad_fn.do ks_multicg_offset.fo ks_multicg_offset.do ks_boundary.fo ks_boundary.do ks_dslash_complete_specialization.fo ks_dslash_complete_specialization.do

OBJECTS+= gauge_force_imp.fo gauge_force_imp.do gf_boundary.fo gf_boundary.do gauge_force_imp_complete_specialization.fo gauge_force_imp_complete_specialization.do

OBJECTS+= ff_boundary.do fermion_force_complete_specialization.do

OBJECTS+= layout.do fermion_force.do qphix_su3_algebra.do #layout.fo fermion_force.fo qphix_su3_algebra.fo

libqphix: ks_dslash_complete_specialization.h ${LIBQPHIXMILC_SRC} $(MAKEFILE) $(CUSTOM_MAKE) ${OBJECTS} fermion_force_complete_specialization.h
	#$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=1 -DVECLEN=${VECLENF} ks_long_dslash.cpp -c  $(CMLLDFLAGS)
	#$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) ks_blas_utils_c.cpp -c  $(CMLLDFLAGS)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=1 -DVECLEN=${VECLENF} qphix_f3.cpp -c $(CMLLDFLAGS)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -DQPHIX_PrecisionInt=2 -DVECLEN=${VECLEND} qphix_d3.cpp -c $(CMLLDFLAGS)
	#$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) qphix_misc.cpp -c $(CMLLDFLAGS)
	#$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) ks_d_congrad_fn_two_src.cpp -c  $(CMLLDFLAGS)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) hrtimer.cpp -c  $(CMLLDFLAGS)
	ar rvs ${QPHIXLIB} ks_long_dslash.fo ks_long_dslash.do ks_blas_utils_c.fo ks_blas_utils_c.do ks_d_congrad_fn.fo ks_d_congrad_fn.do ks_multicg_offset.fo ks_multicg_offset.do qphix_f3.o qphix_d3.o hrtimer.o qphix_misc.fo qphix_misc.do ks_boundary.fo ks_boundary.do ks_dslash_complete_specialization.fo ks_dslash_complete_specialization.do gauge_force_imp.fo gauge_force_imp.do gf_boundary.fo gf_boundary.do ff_boundary.do gauge_force_imp_complete_specialization.fo gauge_force_imp_complete_specialization.do  fermion_force_complete_specialization.do layout.do fermion_force.do qphix_su3_algebra.do  #layout.fo fermion_force.fo qphix_su3_algebra.fo

.PHONY: clean
clean: 
	rm -rf *.o *.fo *.do

.PHONY: install
install:
	mkdir -p $(PREFIX)/include
	mkdir -p $(PREFIX)/lib
	cp $(BASE_DIR)/include/qphix_ff_interface.h $(PREFIX)/include
	cp $(QPHIXLIB) $(PREFIX)/lib
	if test -e $(PREFIX)/lib/libqphixmilc$(SUFFIX).a; then rm $(PREFIX)/lib/libqphixmilc$(SUFFIX).a; fi
	ln -s $(LIB_FILE) $(PREFIX)/lib/libqphixmilc$(SUFFIX).a

.PHONY: uninstall
uninstall:
	rm -rf $(PREFIX)/include
	rm -rf $(PREFIX)/lib
	rmdir $(PREFIX)

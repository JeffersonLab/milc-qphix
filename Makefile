TARGET=ks_long_dslash
#DEBUG=1
mode=mic 
#mode=avx

#mode=sim # Compile for Simulator
#mode=knlsim # Compile for KNL Simulator
#mode=bsd # Compile for FBSD Stack

#CODEGEN_PATH=./

mode:=$(strip $(mode))

CONFFILE=customMake.$(mode)
include $(CONFFILE)

ifeq ($(ENABLE_MPI),1)
ifneq ($(ASSUME_MULTINODE),1)
CXX = mpiicpc
else
CXX = icpc
endif
else
CXX = icpc
endif

ifeq ($(DEBUG),1)
OPTFLAGS = -O0
else
OPTFLAGS = -O3
endif

CXXFLAGS = -Wall -Wextra -g $(OPTFLAGS) -openmp -fargument-noalias-global -restrict $(PREFFLAGS) -I$(SEPHOME)/include -I../qphix-codegen

ifneq ($(knl),1)
SIM_ICCDIR = /p/ctg/arl065/LrbWorkloads/compilers/20110420-lrb2-linux/Compiler
SEE_DIR = /nfs/sc/proj/ctg/arl001/larry/SEE_4.1
SIMFLAGS = -xBETA_R2 -mGLOB_default_function_attrs="use_gather_scatter_hint=off"
SIMLDPATH = -L$(SIM_ICCDIR)/lib/mic2
else
#SIM_ICCDIR = /nfs/sc/proj/ctg/arl062/vkumar7/LrbWorkloads/compilers/latest-avx3-linux
#SEE_DIR = /nfs/sc/proj/ctg/arl062/vkumar7/LrbWorkloads/see/SEE-1.5.2-SSE
SIM_ICCDIR = /nfs/site/disks/lrb0137/lrb2/work/arch01/ddkalamk/20121211-avx3-linux
SEE_DIR = /nfs/site/disks/lrb0137/lrb2/work/arch01/ddkalamk/SEE-1.5.2-SSE
CXXFLAGS = -g $(OPTFLAGS) -fargument-noalias-global $(PREFFLAGS) -I$(SEPHOME)/include
SIMFLAGS = -xBETA_R3 
SIMLDPATH = -L$(SIM_ICCDIR)/lib/mic
DEFS += -DKNLSIM
endif

SIMCXX = $(SIM_ICCDIR)/bin/icpc
SIMINC = -I$(SEE_DIR)/include -I$(SEE_DIR)/include/remap 
SIMLDFLAGS = -N -Ttext 0 -static -L$(SEE_DIR)/lib64 $(SEE_DIR)/lib64/crt1_fish.o $(SEE_DIR)/lib64/crtend_fish.o -L$(SIM_ICCDIR)/lib
SIMLDFLAGS1 = $(SIMLDPATH) -lseekrn -lstdc++ -lsupc++ -lgcc_eh -lmcrtsim -lseekrn -lsvml -lm 

#ifeq ($(mode),mic)
#ifeq ($(PRECISION),1)
#override VECLEN=16
#else
#override VECLEN=8
#endif
#endif

ifeq ($(mode),avx)
ifeq ($(PRECISION),1)
override VECLEN=8
else
override VECLEN=4
endif
endif

ifeq ($(mode),sse)
ifeq ($(PRECISION),1)
override VECLEN=4
else
override VECLEN=2
endif
MINCT=2
MPSSFLAGS = -msse3
endif

ifeq ($(mode),scalar)
override VECLEN=1
override SOALEN=1
endif

ifeq ($(mode),sim)
ifeq ($(PRECISION),1)
override VECLEN=16
else
override VECLEN=8
endif
endif

ifeq ($(mode),mic)
ifneq ($(AVX512),1)
MPSSFLAGS = -mmic -mGLOB_default_function_attrs="use_gather_scatter_hint=off"
else
MPSSFLAGS = -xMIC-AVX512
yesnolist += AVX512
endif
endif

ifeq ($(mode),avx)
ifneq ($(AVX2),1)
MPSSFLAGS = -mavx
else
MPSSFLAGS = -march=core-avx2
endif
MINCT=2
yesnolist += AVX2
endif

ifeq ($(mode),bsd)
CXX = /nfs/sc/proj/ctg/arl065/LrbWorkloads0/compilers/latest-lrb2-linux-bsd/bin/icpc
MPSSFLAGS = -xBETA_R2 -mGLOB_default_function_attrs="use_gather_scatter_hint=off"
DEFS += -DUSE_OLD_INTRIN
TARGET = qcd qcd.s 
ENABLE_SEP_SAMPLING = 0
endif

ifeq ($(mode),sim)
DEFS += -DSIMULATOR
ifneq ($(knl),1)
DEFS += -DUSE_OLD_INTRIN
endif
ENABLE_STREAMING_STORES=0
NREPEAT=8 
BENCH_REPEAT=1
TARGET = qcd.sim
endif

ifeq ($(mode),knlsim)
SIM_ICCDIR = /nfs/sc/proj/ctg/arl062/vkumar7/LrbWorkloads/compilers/latest-avx3-linux
DEFS += -DSIMULATOR
DEFS += -DKNLSIM -DSDE
#DEFS += -DUSE_OLD_INTRIN
NCORES=1
SY=1
SZ=1
NREPEAT=2 
BENCH_REPEAT=1
ENABLE_STREAMING_STORES=0
TARGET = qcd.knlsim
endif

ifeq ($(USE_CML_TG),1)
CMLCXXFLAGS= -I../../synk/trunk -I../../cml-proxy/
CMLLDFLAGS= -lcmlthreading -lsynk -L../../synk/trunk -L../../cml-proxy/
endif

export

yesnolist += USE_PACKED_GAUGES
yesnolist += ENABLE_STREAMING_STORES
yesnolist += ENABLE_LOW_PRECISION
yesnolist += ENABLE_SEP_SAMPLING
yesnolist += ENABLE_BARRIER
yesnolist += COMPRESSED_12
yesnolist += YZ_TILING
yesnolist += CLAMP_OFFCORE
yesnolist += CLAMP_IN_L1
yesnolist += DEBUG_OVRHD
yesnolist += NO_COMPUTE
yesnolist += USE_CML_TG
yesnolist += USE_WAITANY
yesnolist += ENABLE_MPI
yesnolist += ASSUME_MULTINODE
yesnolist += REVERSE_MAP
yesnolist += TIME_CG

deflist += NX
deflist += NY
deflist += NZ
deflist += NT
deflist += NCORES
deflist += MINCT
deflist += SY
deflist += SZ
deflist += BY
deflist += BZ
deflist += NREPEAT
deflist += BENCH_REPEAT
deflist += BARRIER_FREQ_T
deflist += SOALEN
#deflist += VECLEN
#deflist += PRECISION
deflist += XY_PAD
deflist += XYZ_PAD
deflist += CODEGEN_PATH


DEFS += $(strip $(foreach var, $(yesnolist), $(if $(filter 1, $($(var))), -D$(var))))
DEFS += $(strip $(foreach var, $(deflist), $(if $($(var)), -D$(var)=$($(var)))))

LIBQPHIXMILC_SRC= ks_long_dslash.cpp ks_d_congrad_fn.cpp  ks_d_congrad_fn_two_src.cpp ks_multicg_offset.cpp ks_blas_utils_c.cpp hrtimer.cpp                 

all: $(TARGET)

qcd.knlsim: dslash_mic_complete_specialization.h qcd.cpp Makefile customMake.$(mode)
	$(SIMCXX) $(OPTFLAGS) -xBETA_R3 -g -fasm-blocks $(DEFS) qcd.cpp -o qcd.knlsim

qcd.sim: qcd.o
	ld $(SIMLDFLAGS) qcd.o $(SIMLDFLAGS1) -o qcd.sim -liomp5

qcd.o: dslash_mic_complete_specialization.h qcd.cpp Makefile customMake.$(mode)
	$(SIMCXX) $(SIMFLAGS) $(CXXFLAGS) $(DEFS) $(SIMINC) -c qcd.cpp

qcd: dslash_mic_complete_specialization.h qcd.cpp Makefile customMake.$(mode)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) qcd.cpp -o ./qcd $(CMLLDFLAGS)

ks_qcd: ks_dslash_complete_specialization.h ks_qcd.cpp Makefile customMake.$(mode)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) ks_qcd.cpp -o ./ks_qcd $(CMLLDFLAGS)

ks_qcd_imp: ks_dslash_complete_specialization.h ks_qcd_imp.cpp Makefile customMake.$(mode)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) ks_qcd_imp.cpp -o ./ks_qcd_imp $(CMLLDFLAGS)

ks_long_dslash: ks_dslash_complete_specialization.h ${LIBQPHIXMILC_SRC} Makefile customMake.$(mode)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) ks_long_dslash.cpp -c  $(CMLLDFLAGS)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) ks_blas_utils_c.cpp -c  $(CMLLDFLAGS)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) ks_d_congrad_fn.cpp -c  $(CMLLDFLAGS)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) ks_d_congrad_fn_two_src.cpp -c  $(CMLLDFLAGS)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) ks_multicg_offset.cpp -c $(CMLLDFLAGS)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) hrtimer.cpp -c  $(CMLLDFLAGS)
	ar rvs libqphixmilc.a ks_long_dslash.o ks_blas_utils_c.o ks_d_congrad_fn.o ks_multicg_offset.o ks_d_congrad_fn_two_src.o hrtimer.o


qcd.s: dslash_mic_complete_specialization.h qcd.cpp Makefile customMake.$(mode)
	$(CXX) $(CXXFLAGS) $(CMLCXXFLAGS) $(MPSSFLAGS) $(DEFS) -S ./qcd.cpp

run: ./qcd

ifeq ($(mode),mic)
	/opt/intel/mic/coi/tools/micnativeloadex/release/micnativeloadex  ./qcd -d 0 -e "KMP_AFFINITY=compact,granularity=thread"
endif
ifeq ($(mode),avx)
	KMP_AFFINITY=compact ./qcd
endif
ifeq ($(mode),sse)
	KMP_AFFINITY=compact ./qcd
endif
ifeq ($(mode),scalar)
ifeq ($(mic),1)
	/opt/intel/mic/coi/tools/micnativeloadex/release/micnativeloadex  ./qcd -d 0 -e "KMP_AFFINITY=compact,granularity=thread"
else
	KMP_AFFINITY=compact ./qcd
endif
endif

customMake.sim:
	ln -s customMake.mic customMake.sim

copy: ./qcd
	sudo micput mic0 qcd /tmp

cml: synk
	$(MAKE) -C ../../cml-proxy -f Makefile.threading  clean
	$(MAKE) -C ../../cml-proxy -f Makefile.threading

synk:
	$(MAKE) -C ../../synk/trunk clean
	$(MAKE) -C ../../synk/trunk 

clean: 
	rm -rf *.o *.fo *.do ./qcd *.s qcd.sim ./ks_dslash_tester 

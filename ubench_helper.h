#ifndef _UBENCH_HELPER_H_
#define _UBENCH_HELPER_H_

#if (PRECISION == 1) && defined(ENABLE_LOW_PRECISION)
float cvtHalf2Float(half val);
half cvtFloat2Half(float val);
#define DownConv(x) cvtFloat2Half(x)
#define UpConv(x) cvtHalf2Float(x)
#else
#define DownConv(x) (x)
#define UpConv(x) (x)
#endif

#if (PRECISION == 1) && defined(ENABLE_LOW_PRECISION)
#if VECLEN == 16
float cvtHalf2Float(half val) {
	float ret;
	__m512 tmp;
	tmp = _mm512_mask_extloadunpacklo_ps(_mm512_undefined_ps(), 0x1, &val, _MM_UPCONV_PS_FLOAT16, _MM_HINT_NONE);	
	_mm512_mask_packstorelo_ps(&ret, 0x1, tmp);
	return ret;
}
half cvtFloat2Half(float val) {
	half ret;
	__m512 tmp;
	tmp = _mm512_mask_extloadunpacklo_ps(_mm512_undefined_ps(), 0x1, &val, _MM_UPCONV_PS_NONE, _MM_HINT_NONE);	
	_mm512_mask_extpackstorelo_ps(&ret, 0x1, tmp, _MM_DOWNCONV_PS_FLOAT16, _MM_HINT_NONE);
	return ret;
}
#elif VECLEN == 8
float cvtHalf2Float(half val) {
	float ret;
	_mm_store_ss(&ret, _mm_cvtph_ps(_mm_insert_epi16(_mm_setzero_si128(), val, 0)));
	return ret;
}
half cvtFloat2Half(float val) {
	half ret;
	ret =_mm_extract_epi16( _mm_cvtps_ph(_mm_load_ss(&val), _MM_FROUND_CUR_DIRECTION), 0);
	return ret;
}
#endif
#endif

void printHelp() 
{ 
	my_printf("./qcd -x Lx -y Ly -z Lz -t Lt -i iters -by BY -bz BZ -c NCores -sy SY -sz SZ  -pxy xy_pad -pxyz xyz_pad -minct MinCt");
#ifdef ENABLE_MPI
	my_printf(" -geom Px Py Pz Pt");
#ifdef USE_CML_TG
	my_printf(" -parsends nParallelSends -parrecvs nParallelRecvs");
#endif
#endif
	my_printf("\n");
	my_printf("   Lx is the lattice size in X (default %d)\n", NX);
	my_printf("   Ly is the lattice size in Y (default %d)\n", NY);
	my_printf("   Lz is the lattice size in Z (default %d)\n", NZ);
	my_printf("   iters is the number of iterations (default %d)\n", NREPEAT);
	my_printf("   BY is the block size in Y (default %d)\n", BY);
	my_printf("   BZ is the block size in Z (default %d)\n", BZ);
	my_printf("   NCores is the no of cores (default %d)\n", NCORES);
	my_printf("   Sy is the no of simt threads in y (default %d)\n", SY);
	my_printf("   Sz is the no of simt threads in z (default %d)\n", SZ);
	my_printf("   xy_pad is the extra pad in the XY plane (default %d)\n", XY_PAD);
	my_printf("   xyz_pad is the extra pad in the XYZ plane (default %d)\n", XYZ_PAD);
	my_printf("   MinCt is the MinCt parameter in the blocking scheme (default %d)\n", MINCT);
#ifdef ENABLE_MPI
	my_printf("   Px Py Pz Pt define a 4D grid of MPI tasks\n");
	my_printf("   -printall prints communication latency stats for all the ranks\n");
#ifdef USE_CML_TG
	my_printf("   nParallelSends - No. of parallel send groups (1 to 8)\n");
	my_printf("   nParallelRecvs - No. of parallel recv groups (1 or 2)\n");
#endif
#endif
}

void processArgs(int argc, char *argv[]) 
{
	int i=1;
	if (argc == 1) {
		printHelp();
	}

	while( i < argc)  {
		if( strcmp(argv[i], "-x") == 0 ) {
			Gx=atoi(argv[i+1]);
			i+=2;
		}
		else if ( strcmp(argv[i], "-y") == 0 ) {
			Gy=atoi(argv[i+1]);
			i+=2;
		}
		else if ( strcmp(argv[i], "-z") == 0) {
			Gz=atoi(argv[i+1]);
			i+=2;
		}
		else if ( strcmp(argv[i], "-t") == 0) {
			Gt=atoi(argv[i+1]);
			i+=2;
		}
		else if ( strcmp(argv[i], "-i") == 0) {
			iters=atoi(argv[i+1]);
			i+=2;
		}
		else if ( strcmp(argv[i], "-bx") == 0 ) {
			Bx=atoi(argv[i+1]);
			i+=2;
		}
		else if ( strcmp(argv[i], "-by") == 0 ) {
			By=atoi(argv[i+1]);
			i+=2;
		}
		else if (strcmp(argv[i], "-bz") == 0 ) {
			Bz=atoi(argv[i+1]);
			i+=2;
		}
		else if ( strcmp(argv[i], "-c") == 0 ) {
			NCores=atoi(argv[i+1]);
			i+=2;
		}
		else if ( strcmp(argv[i], "-sy") == 0 ) {
			Sy=atoi(argv[i+1]);
			i+=2;
		}
		else if (strcmp(argv[i], "-sz") == 0 ) {
			Sz=atoi(argv[i+1]);
			i+=2;
		}
		else if ( strcmp(argv[i], "-pxy") == 0 ) {
			xy_pad=atoi(argv[i+1]);
			i+=2;
		}
		else if (strcmp(argv[i], "-pxyz") == 0 ) {
			xyz_pad=atoi(argv[i+1]);
			i+=2;
		}
		else if (strcmp(argv[i], "-minct") == 0 ) {
			minCt=atoi(argv[i+1]);
			i+=2;
		}
		else if (strcmp(argv[i], "-dump") == 0 ) {
			dump =true;
			i++;
		}
		else if (strcmp(argv[i], "-printall") == 0 ) {
			printAllRanks =true;
			i++;
		}
		else if (strcmp(argv[i], "-wait") == 0 ) {
			dwait =true;
			waitrank=atoi(argv[i+1]);
			i+=2;
		}
#ifdef ENABLE_MPI
		else if (strcmp(argv[i], "-geom") == 0 ) {
			geometry[0] = atoi(argv[i+1]);
			geometry[1] = atoi(argv[i+2]);
			geometry[2] = atoi(argv[i+3]);
			geometry[3] = atoi(argv[i+4]);
			i+=5;
		}
#ifdef USE_CML_TG
		else if (strcmp(argv[i], "-parsends") == 0 ) {
			nParallelSends=atoi(argv[i+1]);
			i+=2;
		}
		else if (strcmp(argv[i], "-parrecvs") == 0 ) {
			nParallelRecvs=atoi(argv[i+1]);
			i+=2;
		}
#endif
#endif
		else {
			printf("Ignoring arg[%d] = %s\n", i, argv[i]); 
			i++;
		}
	}

	if( NCores < 0 ) { printHelp(); exit(1); }
	if( Sy < 0 ) { printHelp(); exit(1); }
	if( Sz < 0 ) { printHelp(); exit(1); }

	// Ct does not have to divide t, we can pick that up.
	if( By < 0 ) { By = Gy; }
	if( Bz < 0 ) { Bz = Gz; }

	n_threads_per_core = Sy*Sz;
	if(minCt > NCores) minCt = NCores;

}

void checkParams()
{
	MYASSERT(Nxh % 2 == 0);
	MYASSERT(Ny % 4 == 0);
	MYASSERT(Nz % 4 == 0);
	MYASSERT(Nt % 4 == 0);

	Vx = (VECLEN > 1 ? Nx/2 : Nx);
	Vxh = (VECLEN > 1 ? Nxh/2 : Nxh);
	Vy = (VECLEN > 2 ? Ny/2 : Ny);
	Vz = (VECLEN > 4 ? Nz/2 : Nz);
	Vt = (VECLEN > 8 ? Nt/2 : Nt);
	if(Bx == 0) Bx = Vxh;
	//MYASSERT(nvecs % Bx == 0);
	MYASSERT(Vy % By == 0);
	MYASSERT(Vz % Bz == 0);
	Pxy = (Vxh*Vy+xy_pad);
	Pxyz = (Pxy*Vz+xyz_pad);
}

void debug_wait()
{
	if(dwait && (waitrank == -1 || waitrank == myRank)) {
		printf("Rank %d (pid=%d) Waiting on dummy\n", myRank, getpid());
		volatile int dummy = 0;
		while(dummy == 0) ;
	}
}

void init_spinor(Spinor *s_in)
{
	int nGX = Nxh/Vxh;
	int nGY = Ny/Vy;
	int nGZ = Nz/Vz;
	int nGT = Nt/Vt;

#ifdef _OPENMP
#pragma omp parallel for collapse(4)
#endif
	for(int t = 0; t < Vt; t++) {
		for(int z = 0; z < Vz; z++) {
			for(int y = 0; y < Vy; y++) {
				for(int x = 0; x < Vxh; x++) {
					int ind = t*Pxyz+z*Pxy+y*Vxh+x;
					for(int sp = 0; sp < 4; sp++) {
						for(int c = 0; c < 3; c++) {
							for(int gt = 0; gt < nGT; gt++) {
								for(int gz = 0; gz < nGZ; gz++) {
									for(int gy = 0; gy < nGY; gy++) {
										for(int gx = 0; gx < nGX ; gx++) {
											int v = ((gt*nGZ+gz)*nGY+gy)*nGX+gx;
											s_in[ind][c][sp][0][v] = DownConv((gt*Vt+t+Lst)*0.1+(gz*Vz+z+Lsz)*0.2-(gy*Vy+y+Lsy)*0.3+(gx*Vxh+x+Lsxh)*0.4-sp*0.1+c*0.2-0.3);
											s_in[ind][c][sp][1][v] = DownConv(-((gt*Vt+t+Lst)*0.1-(gz*Vz+z+Lsz)*0.2+(gy*Vy+y+Lsy)*0.3-(gx*Vxh+x+Lsxh)*0.4+sp*0.1-c*0.2+0.4));
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

void init_spinor(KS * s_in)
{
	int nGX = Nxh/Vxh;
	int nGY = Ny/Vy;
	int nGZ = Nz/Vz;
	int nGT = Nt/Vt;

#ifdef _OPENMP
#pragma omp parallel for collapse(4)
#endif
	for(int t = 0; t < Vt; t++) {
		for(int z = 0; z < Vz; z++) {
			for(int y = 0; y < Vy; y++) {
				for(int x = 0; x < Vxh; x++) {
					int ind = t*Pxyz+z*Pxy+y*Vxh+x;
					for(int c = 0; c < 3; c++) {
						for(int gt = 0; gt < nGT; gt++) {
							for(int gz = 0; gz < nGZ; gz++) {
								for(int gy = 0; gy < nGY; gy++) {
									for(int gx = 0; gx < nGX ; gx++) {
										int v = ((gt*nGZ+gz)*nGY+gy)*nGX+gx;
										s_in[ind][c][0][v] = DownConv((gt*Vt+t+Lst)*0.1+(gz*Vz+z+Lsz)*0.2-(gy*Vy+y+Lsy)*0.3+(gx*Vxh+x+Lsxh)*0.4+c*0.2-0.3);
										s_in[ind][c][1][v] = DownConv(-((gt*Vt+t+Lst)*0.1-(gz*Vz+z+Lsz)*0.2+(gy*Vy+y+Lsy)*0.3-(gx*Vxh+x+Lsxh)*0.4-c*0.2+0.4));
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

void init_gauge(Gauge * gfl)
{
	int nGX = Nxh/Vxh;
	int nGY = Ny/Vy;
	int nGZ = Nz/Vz;
	int nGT = Nt/Vt;

#ifdef _OPENMP
#pragma omp parallel for collapse(4)
#endif
	for(int t = 0; t < Vt; t++) {
		for(int z = 0; z < Vz; z++) {
			for(int y = 0; y < Vy; y++) {
				for(int x = 0; x < Vxh; x++) {
					int ind = t*Pxyz+z*Pxy+y*Vxh+x;
					for(int d = 0; d < 8; d++) {
						for(int c = 0; c < (COMPRESSED_GAUGES ? 2 : 3); c++) {
							for(int c2 = 0; c2 < 3; c2++) {
								for(int gt = 0; gt < nGT; gt++) {
									for(int gz = 0; gz < nGZ; gz++) {
										for(int gy = 0; gy < nGY; gy++) {
											for(int gx = 0; gx < nGX ; gx++) {
												int v = ((gt*nGZ+gz)*nGY+gy)*nGX+gx;
												gfl[ind][d][c][c2][0][v] = DownConv(4.0e-3*((gt*Vt+t+Lst)*0.1-(gz*Vz+z+Lsz)*0.2+(gy*Vy+y+Lsy)*0.3-(gx*Vxh+x+Lsxh)*0.4+c*0.1-c2*0.2+0.3*d-0.5));
												gfl[ind][d][c][c2][1][v] = DownConv(-4.0e-3*((gt*Vt+t+Lst)*0.1+(gz*Vz+z+Lsz)*0.2-(gy*Vy+y+Lsy)*0.3+(gx*Vxh+x+Lsxh)*0.4-c*0.1+c2*0.2-0.3*d+0.4));
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

void init_gauge18(Gauge18 * gfl)
{
	int nGX = Nxh/Vxh;
	int nGY = Ny/Vy;
	int nGZ = Nz/Vz;
	int nGT = Nt/Vt;

#ifdef _OPENMP
#pragma omp parallel for collapse(4)
#endif
	for(int t = 0; t < Vt; t++) {
		for(int z = 0; z < Vz; z++) {
			for(int y = 0; y < Vy; y++) {
				for(int x = 0; x < Vxh; x++) {
					int ind = t*Pxyz+z*Pxy+y*Vxh+x;
					for(int d = 0; d < 8; d++) {
						for(int c = 0; c < 3; c++) {
							for(int c2 = 0; c2 < 3; c2++) {
								for(int gt = 0; gt < nGT; gt++) {
									for(int gz = 0; gz < nGZ; gz++) {
										for(int gy = 0; gy < nGY; gy++) {
											for(int gx = 0; gx < nGX ; gx++) {
												int v = ((gt*nGZ+gz)*nGY+gy)*nGX+gx;
												gfl[ind][d][c][c2][0][v] = DownConv(9.0e-3*((gt*Vt+t+Lst)*0.1-(gz*Vz+z+Lsz)*0.2+(gy*Vy+y+Lsy)*0.3-(gx*Vxh+x+Lsxh)*0.4+c*0.1-c2*0.2+0.3*d-0.5));
												gfl[ind][d][c][c2][1][v] = DownConv(-9.0e-3*((gt*Vt+t+Lst)*0.1+(gz*Vz+z+Lsz)*0.2-(gy*Vy+y+Lsy)*0.3+(gx*Vxh+x+Lsxh)*0.4-c*0.1+c2*0.2-0.3*d+0.4));
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

void set_zero(fptype *buf, int size)
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
#pragma simd
	for(int i=0; i < size; i++)
	{
		buf[i]=0.0;
	}
}

void checksum_spinor(Spinor *si, char *name)
{
	double checksum=0;
	int nGX = Nxh/Vxh;
	int nGY = Ny/Vy;
	int nGZ = Nz/Vz;
	int nGT = Nt/Vt;

#ifdef _OPENMP
#pragma omp parallel for collapse(4) reduction(+:checksum)
#endif
	for(int t = 0; t < Vt; t++) {
		for(int z = 0; z < Vz; z++) {
			for(int y = 0; y < Vy; y++) {
				for(int x = 0; x < Vxh; x++) {
					int ind = t*Pxyz+z*Pxy+y*Vxh+x;
					for(int sp = 0; sp < 4; sp++) {
						for(int c = 0; c < 3; c++) {
							for(int gt = 0; gt < nGT; gt++) {
								for(int gz = 0; gz < nGZ; gz++) {
									for(int gy = 0; gy < nGY; gy++) {
										for(int gx = 0; gx < nGX ; gx++) {
											int v = ((gt*nGZ+gz)*nGY+gy)*nGX+gx;
											checksum += UpConv(si[ind][c][sp][0][v]) * ((gt*Vt+t+Lst)*0.1+(gz*Vz+z+Lsz)*0.2-(gy*Vy+y+Lsy)*0.3+(gx*Vxh+x+Lsxh)*0.4-sp*0.1+c*0.2-0.3);
											checksum += UpConv(si[ind][c][sp][1][v]) * ((gt*Vt+t+Lst)*0.1-(gz*Vz+z+Lsz)*0.2+(gy*Vy+y+Lsy)*0.3-(gx*Vxh+x+Lsxh)*0.4+sp*0.1-c*0.2+0.4);
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	//printf("Rank %d: KS Local checksum=%lg\n", myRank, checksum);
	checksum = mpi_allreduce(checksum);
	my_printf("Spinor (%s) checksum=%lg\n", name, checksum);
}

void dump_spinor(Spinor *si)
{
	int nGX = Nxh/Vxh;
	int nGY = Ny/Vy;
	int nGZ = Nz/Vz;
	int nGT = Nt/Vt;

	for(int t = 0; t < Vt; t++) {
		for(int z = 0; z < Vz; z++) {
			for(int y = 0; y < Vy; y++) {
				for(int x = 0; x < Vxh; x++) {
					int ind = t*Pxyz+z*Pxy+y*Vxh+x;
					for(int sp = 0; sp < 4; sp++) {
						for(int c = 0; c < 3; c++) {
							for(int gt = 0; gt < nGT; gt++) {
								for(int gz = 0; gz < nGZ; gz++) {
									for(int gy = 0; gy < nGY; gy++) {
										for(int gx = 0; gx < nGX ; gx++) {
											int v = ((gt*nGZ+gz)*nGY+gy)*nGX+gx;
										  int gind = (gx*Vxh+x+Lsxh) + (Gx/2) * ((gy*Vy+y+Lsy)+Gy*((gz*Vz+z+Lsz)+Gz*((gt*Vt+t+Lst))));
											printf("%05d %d %d %d %d %d %d %g %g -- BLOCKS\n", gind, (gt*Vt+t+Lst), (gz*Vz+z+Lsz), (gy*Vy+y+Lsy), (gx*Vxh+x+Lsxh), sp, c, UpConv(si[ind][c][sp][0][v]), UpConv(si[ind][c][sp][1][v]));
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

void checksum_spinor(KS *si, char* name)
{
	double checksum=0;
	int nGX = Nxh/Vxh;
	int nGY = Ny/Vy;
	int nGZ = Nz/Vz;
	int nGT = Nt/Vt;

#ifdef _OPENMP
#pragma omp parallel for collapse(4) reduction(+:checksum)
#endif
	for(int t = 0; t < Vt; t++) {
		for(int z = 0; z < Vz; z++) {
			for(int y = 0; y < Vy; y++) {
				for(int x = 0; x < Vxh; x++) {
					int ind = t*Pxyz+z*Pxy+y*Vxh+x;
					for(int c = 0; c < 3; c++) {
						for(int gt = 0; gt < nGT; gt++) {
							for(int gz = 0; gz < nGZ; gz++) {
								for(int gy = 0; gy < nGY; gy++) {
									for(int gx = 0; gx < nGX ; gx++) {
										int v = ((gt*nGZ+gz)*nGY+gy)*nGX+gx;
										checksum += UpConv(si[ind][c][0][v]) * ((gt*Vt+t+Lst)*0.1+(gz*Vz+z+Lsz)*0.2-(gy*Vy+y+Lsy)*0.3+(gx*Vxh+x+Lsxh)*0.4+c*0.2-0.3);
										checksum += UpConv(si[ind][c][1][v]) * ((gt*Vt+t+Lst)*0.1-(gz*Vz+z+Lsz)*0.2+(gy*Vy+y+Lsy)*0.3-(gx*Vxh+x+Lsxh)*0.4-c*0.2+0.4);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	//printf("Rank %d: KS Local checksum=%lg\n", myRank, checksum);
	checksum = mpi_allreduce(checksum);
	my_printf("KS (%s) checksum=%lg\n", name, checksum);
}

void dump_spinor(KS *si)
{
	int nGX = Nxh/Vxh;
	int nGY = Ny/Vy;
	int nGZ = Nz/Vz;
	int nGT = Nt/Vt;

	for(int t = 0; t < Vt; t++) {
		for(int z = 0; z < Vz; z++) {
			for(int y = 0; y < Vy; y++) {
				for(int x = 0; x < Vxh; x++) {
					int ind = t*Pxyz+z*Pxy+y*Vxh+x;
					for(int c = 0; c < 3; c++) {
						for(int gt = 0; gt < nGT; gt++) {
							for(int gz = 0; gz < nGZ; gz++) {
								for(int gy = 0; gy < nGY; gy++) {
									for(int gx = 0; gx < nGX ; gx++) {
										int v = ((gt*nGZ+gz)*nGY+gy)*nGX+gx;
										int gind = (gx*Vxh+x+Lsxh) + (Gx/2) * ((gy*Vy+y+Lsy)+Gy*((gz*Vz+z+Lsz)+Gz*((gt*Vt+t+Lst))));
										printf("%05d %d %d %d %d %d %g %g -- BLOCKS\n", gind, (gt*Vt+t+Lst), (gz*Vz+z+Lsz), (gy*Vy+y+Lsy), (gx*Vxh+x+Lsxh), c, UpConv(si[ind][c][0][v]), UpConv(si[ind][c][1][v]));
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

void checksum_gauge(Gauge *gfl, char *name)
{
	double checksum=0;
	int nGX = Nxh/Vxh;
	int nGY = Ny/Vy;
	int nGZ = Nz/Vz;
	int nGT = Nt/Vt;

#ifdef _OPENMP
#pragma omp parallel for collapse(4) reduction(+:checksum)
#endif
	for(int t = 0; t < Vt; t++) {
		for(int z = 0; z < Vz; z++) {
			for(int y = 0; y < Vy; y++) {
				for(int x = 0; x < Vxh; x++) {
					int ind = t*Pxyz+z*Pxy+y*Vxh+x;
					for(int d = 0; d < 8; d++) {
						for(int c = 0; c < (COMPRESSED_GAUGES ? 2 : 3); c++) {
							for(int c2 = 0; c2 < 3; c2++) {
								for(int gt = 0; gt < nGT; gt++) {
									for(int gz = 0; gz < nGZ; gz++) {
										for(int gy = 0; gy < nGY; gy++) {
											for(int gx = 0; gx < nGX ; gx++) {
												int v = ((gt*nGZ+gz)*nGY+gy)*nGX+gx;
												checksum += UpConv(gfl[ind][d][c][c2][0][v]) * 1.0e-3*((gt*Vt+t+Lst)*0.1-(gz*Vz+z+Lsz)*0.2+(gy*Vy+y+Lsy)*0.3-(gx*Vxh+x+Lsxh)*0.4+c*0.1-c2*0.2+0.3*d-0.5);
												checksum += UpConv(gfl[ind][d][c][c2][1][v]) * 1.0e-3*((gt*Vt+t+Lst)*0.1+(gz*Vz+z+Lsz)*0.2-(gy*Vy+y+Lsy)*0.3+(gx*Vxh+x+Lsxh)*0.4-c*0.1+c2*0.2-0.3*d+0.4);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	checksum = mpi_allreduce(checksum);
	my_printf("Gauge (%s) checksum=%lg\n", name, checksum);
}

void checksum_gauge18(Gauge18 *gfl, char *name)
{
	double checksum=0;
	int nGX = Nxh/Vxh;
	int nGY = Ny/Vy;
	int nGZ = Nz/Vz;
	int nGT = Nt/Vt;

#ifdef _OPENMP
#pragma omp parallel for collapse(4) reduction(+:checksum)
#endif
	for(int t = 0; t < Vt; t++) {
		for(int z = 0; z < Vz; z++) {
			for(int y = 0; y < Vy; y++) {
				for(int x = 0; x < Vxh; x++) {
					int ind = t*Pxyz+z*Pxy+y*Vxh+x;
					for(int d = 0; d < 8; d++) {
						for(int c = 0; c < 3; c++) {
							for(int c2 = 0; c2 < 3; c2++) {
								for(int gt = 0; gt < nGT; gt++) {
									for(int gz = 0; gz < nGZ; gz++) {
										for(int gy = 0; gy < nGY; gy++) {
											for(int gx = 0; gx < nGX ; gx++) {
												int v = ((gt*nGZ+gz)*nGY+gy)*nGX+gx;
												checksum += UpConv(gfl[ind][d][c][c2][0][v]) * 1.0e-3*((gt*Vt+t+Lst)*0.1-(gz*Vz+z+Lsz)*0.2+(gy*Vy+y+Lsy)*0.3-(gx*Vxh+x+Lsxh)*0.4+c*0.1-c2*0.2+0.3*d-0.5);
												checksum += UpConv(gfl[ind][d][c][c2][1][v]) * 1.0e-3*((gt*Vt+t+Lst)*0.1+(gz*Vz+z+Lsz)*0.2-(gy*Vy+y+Lsy)*0.3+(gx*Vxh+x+Lsxh)*0.4-c*0.1+c2*0.2-0.3*d+0.4);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	checksum = mpi_allreduce(checksum);
	my_printf("Gauge18 (%s) checksum=%lg\n", name, checksum);
}
#endif // _UBENCH_HELPER_H_

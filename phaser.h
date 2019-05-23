#ifndef _QPHIX_PHASER_H_
#define _QPHIX_PHASER_H_

#define MAX_PHASES 8

class Phaser 
{
	struct CorePhase {
		int Ct;
		int Cyz;
		int startBlock;
		int phRepeats;
	};
	struct BlockPhase {
		int bx;
		int by;
		int bz;
		int bt;
		int nt;
		int mainPhs;
	};

	struct ThreadBlockInfo {
		int group_tid;
		int cid_t;
		int cid_yz;
	};

private:
	CorePhase phases[MAX_PHASES];
	int PHS, repPhases;
	Barrier*** barriers;
	BlockPhase **block_info;
	ThreadBlockInfo **tb_info; 
	int *cphr;
	int lx, ly, lz;

	int Nxh, Ny, Nz, Nt, Bx, By, Bz, NCores, Sy, Sz, minCt, n_threads_per_core;

public:
	struct PhaserState {
		int ph, ph1, cx, cy, cz, ct, Nct;
		int tid, cid, smtid, smtid_y, smtid_z;
	};

	Phaser(int Nxh_, int Ny_, int Nz_, int Nt_, int Bx_, int By_, int Bz_, int NCores_, int Sy_, int Sz_, int minCt_) :
		Nxh(Nxh_), Ny(Ny_), Nz(Nz_), Nt(Nt_), Bx(Bx_), By(By_), Bz(Bz_), NCores(NCores_), Sy(Sy_), Sz(Sz_), minCt(minCt_), n_threads_per_core(Sy*Sz)
	{
		PHS = 0;
		repPhases = 0;

		lx = Nxh / Bx;
		ly = Ny / By;
		lz = Nz / Bz;
		int rem = lx * ly * lz;
		int stblk = 0;
		int nCores = NCores / minCt;
		while(rem > 0) {
			int ctd = nCores / rem;
			int ctu = (nCores + rem - 1) / rem;
			CorePhase &phase = phases[PHS];
			phase.Ct = (ctu <= 4 ? ctu : ctd) * minCt;
			if(phase.Ct > Nt) phase.Ct = Nt;
			phase.Cyz = NCores / phase.Ct;
			if(phase.Cyz > rem) phase.Cyz = rem;
			phase.startBlock = stblk;
			phase.phRepeats = 0;
			do {
				stblk += phase.Cyz;
				rem -= phase.Cyz;
				phase.phRepeats++;
			} while(rem >= phase.Cyz);
			/*printf("Phase %d: Cyz = %d Ct = %d, start = %d, repeat = %d\n"
                   , PHS, phase.Cyz, phase.Ct
                   , phase.startBlock, phase.phRepeats);*/
                   repPhases += phase.phRepeats;
			PHS++;
			if(PHS > MAX_PHASES) {
				printf("PHS run out of space - Please increase MAX_PHASES and recompile\n");
				exit(1);
			}
		}
		barriers = new Barrier**[PHS];
		for(int ph = 0; ph < PHS; ph++) {
			barriers[ph] = new Barrier*[phases[ph].Ct];
			for(int i=0; i < phases[ph].Ct; i++) 
                barriers[ph][i]=new Barrier(phases[ph].Cyz,Sy*Sz);
		}
		block_info = new BlockPhase*[NCores];
		tb_info = new ThreadBlockInfo*[NCores*Sy*Sz];
		cphr = new int[NCores];
	}

	void init(int tid)
	{
		int cid = tid / n_threads_per_core;
		int smtid = tid - n_threads_per_core * cid;
		int phr = 0;
		if(smtid == 0) {
			block_info[cid] = new BlockPhase[repPhases];
			tb_info[cid] = new ThreadBlockInfo[PHS];
		}
		for(int ph = 0; ph < PHS; ph++) {
			CorePhase &phase = phases[ph];
			int nActiveCores = phase.Cyz * phase.Ct;
			if(cid >= nActiveCores) continue;
			int cid_t = cid/phase.Cyz;
			int ngroup = phase.Cyz * n_threads_per_core;
			int group_tid = tid % ngroup;
			if(smtid == 0) {
				ThreadBlockInfo &ti = tb_info[cid][ph];

				ti.cid_t = cid_t;
				ti.cid_yz = cid - cid_t * phase.Cyz;
				//ti.group_tid = group_tid;

				int syz = phase.startBlock + ti.cid_yz;
				for(int r = 0; r < phase.phRepeats; r++) {
					BlockPhase &bi = block_info[cid][phr];
					int tmp1 = syz / lx;
					bi.bx = syz - tmp1 * lx;
					bi.bz = tmp1 / ly;
					bi.by = tmp1 - bi.bz * ly;
					bi.bt = (Nt*ti.cid_t) / phase.Ct;
					bi.nt = (Nt*(ti.cid_t+1)) / phase.Ct - bi.bt;
					//my_printf("CID %02d: TID %03d: bt[%d]=%d, nt[%d]=%d\n", cid, tid, phr, bt[phr], phr, nt[phr]);
					bi.bx *= Bx;
					bi.by *= By;
					bi.bz *= Bz;
					bi.mainPhs = ph;
					syz += phase.Cyz;
					phr++;
				}
			}
			barriers[ph][cid_t]->init(group_tid);
		}
		if(smtid==0) cphr[cid] = phr;
	}

	bool start(PhaserState &ps, const int tid, int &x, int &y, int &z, int &t)
	{
        
		int cid = tid / n_threads_per_core;
		int smtid = tid - n_threads_per_core*cid;
		//printf("DEBUG: Phase.start. tid=%d cid=%d\n", tid, cid);
        	// Compute smt ID Y and Z indices
		int smtid_z = smtid / Sy;
		int smtid_y = smtid - Sy * smtid_z;
		if(cphr[cid] == 0) return false;
		ps.cid = cid;
		ps.tid = tid;
		ps.smtid = smtid;
		ps.smtid_y = smtid_y;
		ps.smtid_z = smtid_z;
		ps.ph = 0;
		BlockPhase &bi = block_info[cid][ps.ph];
        	//printf("DEBUG: Phase.start. bi.bt=%d bi.bz=%d\n", bi.bt , bi.bz);
		ps.ph1 = bi.mainPhs;
		ps.cx = 0;
		ps.cy = smtid_y; //nyg*smtid_y;
		ps.cz = smtid_z;
		ps.ct = 0;
		ps.Nct = bi.nt;

		t = bi.bt;
		z = bi.bz + smtid_z;
		y = bi.by + ps.cy;
		x = bi.bx;
		return true;
	}

	bool next(PhaserState &ps, int &x, int &y, int &z, int &t) {
		bool ret = true;
		ps.cx = ps.cx + 1;
		if(ps.cx == Bx) { 
			ps.cx = 0; 
			ps.cy += Sy; //nyg*Sy;
			if(ps.cy >= By) { 
				ps.cy = ps.smtid_y;//nyg*ps.smtid_y; 
				ps.cz += Sz;
				if(ps.cz >= Bz) { 
					ps.cz = ps.smtid_z; 
					ps.ct++;
					if(ps.ct >= ps.Nct) { 
						ps.ct = 0; 
						ps.ph++;
						ps.Nct = block_info[ps.cid][ps.ph].nt;
						if(ps.ph == cphr[ps.cid]) {
							ps.ph = 0;
							ret = false;
						}
					}
				}
			}
		}
		BlockPhase &bi = block_info[ps.cid][ps.ph];
		t = bi.bt + ps.ct;
		z = bi.bz + ps.cz;
		y = bi.by + ps.cy;
		x = bi.bx + ps.cx;

		if((x < 0) || (x >= Nxh) || (y < 0) || (y >= Ny) || (z < 0) || (z >= Nz) || (t < 0) || (t >= Nt)) {
			printf("YYY Tid = %d, xyzt = %2d %2d %2d %2d, ph=%d, Nct=%d, phr=%d\n", ps.tid, x, y, z, t, ps.ph, ps.Nct, cphr[ps.cid]);
			exit(1);
		}

		return ret;
	}

	void sync(PhaserState &ps) {
		int group_tid = tb_info[ps.cid][ps.ph1].cid_yz * Sy *Sz+ps.smtid;
		barriers[ps.ph1][tb_info[ps.cid][ps.ph1].cid_t]->wait(group_tid);

	}

};

#endif


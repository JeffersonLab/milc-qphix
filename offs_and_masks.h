#ifndef _OFFS_AND_MASKS_H_
#define _OFFS_AND_MASKS_H_

#include <cstddef>
class SpinorOnM1
{
	private:
		__declspec(align(64)) int xbOffs_xodd[2][VECLEN], xbOffs_x0_xodd[2][VECLEN];
		__declspec(align(64)) int xfOffs_xodd[2][VECLEN], xfOffs_xn_xodd[2][VECLEN];
		__declspec(align(64)) int ybOffs_yn0[VECLEN], ybOffs_y0[VECLEN], yfOffs_ynn[VECLEN], yfOffs_yn[VECLEN];
		__declspec(align(64)) int offs[VECLEN];

		unsigned int xbmask_x0_xodd[2];
		unsigned int xfmask_xn_xodd[2];
		unsigned int ybmask_y0;
		unsigned int yfmask_yn;

		const int nvecs, Ny;
		const std::size_t spinor_size;
		const int nyg;

	public:
		SpinorOnM1(int nvecs_, int Ny_, std::size_t ss, bool *local_dir) 
            : nvecs(nvecs_), Ny(Ny_), spinor_size(ss), nyg(VECLEN/SOALEN) 
        {
			xbmask_x0_xodd[0] = xbmask_x0_xodd[1] = -1;
			xfmask_xn_xodd[0] = xfmask_xn_xodd[1] = -1;
			ybmask_y0 = yfmask_yn = -1;

			int spinor_line_in_fptypes = ss/sizeof(fptype);
			for(int y = 0; y < nyg; y++) {
				// Various indexing things
				int ind = y*SOALEN;
				int X = nvecs * y * spinor_line_in_fptypes;
				int y1 = y & 1;
				int y2 = 1 - y1;
				for(int x = 0; x < SOALEN; x++) {
					xbOffs_x0_xodd[y1][ind] = X + x - 1;
					xbOffs_xodd[y1][ind] = X + x - 1;
					if(x == 0) {
						if(local_dir[0]) {
							xbOffs_x0_xodd[y1][ind] -= (spinor_line_in_fptypes - SOALEN - nvecs * spinor_line_in_fptypes);
						}
						else {
							xbOffs_x0_xodd[y1][ind] += SOALEN; // This lane is disabled, just set it within same cache line
							xbmask_x0_xodd[y1] &= ~(1 << ind); // reset a bit in the mask
						}
						xbOffs_xodd[y1][ind]    -= (spinor_line_in_fptypes - SOALEN);
					}
					xfOffs_xodd[y1][ind] = X + x;
					xfOffs_xn_xodd[y1][ind] = X + x;

					xbOffs_x0_xodd[y2][ind] = X + x;
					xbOffs_xodd[y2][ind] = X + x;
					xfOffs_xodd[y2][ind] = X + x + 1;
					xfOffs_xn_xodd[y2][ind] = X + x + 1;
					if(x == SOALEN - 1) {
						xfOffs_xodd[y2][ind] += (spinor_line_in_fptypes - SOALEN);
						if(local_dir[0]) {
							xfOffs_xn_xodd[y2][ind] += (spinor_line_in_fptypes - SOALEN - nvecs * spinor_line_in_fptypes);
						}
						else {
							xfOffs_xn_xodd[y2][ind] -= SOALEN; // This lane is disabled, just set it within same cache line
							xfmask_xn_xodd[y2] &= ~(1 << ind); // reset the ind bit in the mask
						}
					}

					ybOffs_y0[ind] = X - nvecs*spinor_line_in_fptypes + x; // previous y-neighbor site offsets
					if(y == 0) {
						if(local_dir[1]) {
							ybOffs_y0[ind] += Ny*nvecs*spinor_line_in_fptypes;
						} 
						else {
							ybOffs_y0[ind] = X + x; // This lane is disabled, just set it within same cache line
							ybmask_y0 &= ~(1 << ind); // reset the ind bit in the mask 
						}
					}
					ybOffs_yn0[ind] = X - nvecs*spinor_line_in_fptypes + x; // previous y-neighbor site offsets
					yfOffs_yn[ind] = X +  nvecs*spinor_line_in_fptypes + x; // next y-neighbor site offsets
					if(y == nyg - 1) {
						if(local_dir[1]) {
							yfOffs_yn[ind] -= Ny*nvecs*spinor_line_in_fptypes;
						} 
						else {
							yfOffs_yn[ind] = X + x; // This lane is disabled, just set it within same cache line
							yfmask_yn &= ~(1 << ind); // reset the ind bit in the mask 
						}
					}
					yfOffs_ynn[ind] = X +  nvecs*spinor_line_in_fptypes + x; // next y-neighbor site offsets
					offs[ind] = X + x;  // site offsets for z & t neighbors
					ind++;
				}
			}
			if(nyg == 1) {
				if(!local_dir[1]) ybmask_y0 = yfmask_yn = 0;
			}
		}
		int*
        getXbOffs(int x, int xodd) 
        {
            return (x == 0 ? xbOffs_x0_xodd[xodd] : xbOffs_xodd[xodd]);
        }
        
		int* 
        getXfOffs(int x, int xodd) 
        {
            return (x == nvecs-1 ? xfOffs_xn_xodd[xodd] : xfOffs_xodd[xodd]);
        }
        
		int* 
        getYbOffs(int y) 
        {   
            return (y == 0 ? ybOffs_y0 : ybOffs_yn0);
        }
        
		int* 
        getYfOffs(int y) 
        {
            return (y == Ny - nyg ? yfOffs_yn : yfOffs_ynn);
        }
        
		int*
        getOffs() 
        {   
            return offs;
        }

		unsigned int 
        getXbMask(int x, int xodd) 
        {
            return (x == 0 ? xbmask_x0_xodd[xodd] : -1);
        }
        
		unsigned int 
        getXfMask(int x, int xodd) 
        {
            return (x == nvecs-1 ? xfmask_xn_xodd[xodd] : -1);
        }
        
		unsigned int 
        getYbMask(int y) 
        {
            return (y == 0 ? ybmask_y0 : -1);
        }
        
		unsigned int 
        getYfMask(int y) 
        {
            return (y == Ny - nyg ? yfmask_yn : -1);
        }
};


class SpinorOnM3
{
	private:
		__declspec(align(64)) int xbOffs_xodd[2][VECLEN], xbOffs_x0_xodd[2][VECLEN];
		__declspec(align(64)) int xfOffs_xodd[2][VECLEN], xfOffs_xn_xodd[2][VECLEN];
		__declspec(align(64)) int ybOffs_yn0[VECLEN], ybOffs_y0[VECLEN], ybOffs_y2[VECLEN], yfOffs_ynn[VECLEN], yfOffs_yn[VECLEN], yfOffs_yn2[VECLEN];
		__declspec(align(64)) int offs[VECLEN];

		unsigned int xbmask_x0_xodd[2];
		unsigned int xfmask_xn_xodd[2];
		unsigned int ybmask_y0;
		unsigned int ybmask_y2;
		unsigned int yfmask_yn;
		unsigned int yfmask_yn2;

		const int nvecs, Ny;
		const std::size_t spinor_size;
		const int nyg;

	public:
		SpinorOnM3(int nvecs_, int Ny_, std::size_t ss, bool *local_dir) 
            : nvecs(nvecs_), Ny(Ny_), spinor_size(ss), nyg(VECLEN/SOALEN) 
        {
			xbmask_x0_xodd[0] = xbmask_x0_xodd[1] = -1;
			xfmask_xn_xodd[0] = xfmask_xn_xodd[1] = -1;
			ybmask_y0 = yfmask_yn = -1;
			ybmask_y2 = yfmask_yn2 = -1;

			int spinor_line_in_fptypes = ss/sizeof(fptype);
			for(int y = 0; y < nyg; y++) {
				// Various indexing things
				int ind = y*SOALEN;
				int X = nvecs * y * spinor_line_in_fptypes;
				int y1 = y & 1;
				int y2 = 1 - y1;
				for(int x = 0; x < SOALEN; x++) {
					xbOffs_x0_xodd[y1][ind] = X + x - 2;
					xbOffs_xodd[y1][ind] = X + x - 2;
					if(x <= 1) {
						if(local_dir[0]) {
							xbOffs_x0_xodd[y1][ind] -= (spinor_line_in_fptypes - SOALEN - nvecs * spinor_line_in_fptypes);
						}
						else {
							xbOffs_x0_xodd[y1][ind] += SOALEN; // This lane is disabled, just set it within same cache line
							xbmask_x0_xodd[y1] &= ~(1 << ind); // reset a bit in the mask
						}
						xbOffs_xodd[y1][ind]    -= (spinor_line_in_fptypes - SOALEN);
					}
					xfOffs_xodd[y1][ind] = X + x + 1;
					xfOffs_xn_xodd[y1][ind] = X + x + 1;

					if(x == SOALEN - 1) {
						if(local_dir[0]) {
							xfOffs_xn_xodd[y1][ind] += (spinor_line_in_fptypes - SOALEN - nvecs * spinor_line_in_fptypes);
						}
						else {
							xfOffs_xn_xodd[y1][ind] -= SOALEN; // This lane is disabled, just set it within same cache line
							xfmask_xn_xodd[y1] &= ~(1 << ind); // reset a bit in the mask
						}
						xfOffs_xodd[y1][ind]    += (spinor_line_in_fptypes - SOALEN);
					}

					xbOffs_x0_xodd[y2][ind] = X + x - 1;
					xbOffs_xodd[y2][ind] = X + x - 1;
					if(x == 0) {
						if(local_dir[0]) {
							xbOffs_x0_xodd[y2][ind] -= (spinor_line_in_fptypes - SOALEN - nvecs * spinor_line_in_fptypes);
						}
						else {
							xbOffs_x0_xodd[y2][ind] += SOALEN; // This lane is disabled, just set it within same cache line
							xbmask_x0_xodd[y2] &= ~(1 << ind); // reset a bit in the mask
						}
						xbOffs_xodd[y2][ind]    -= (spinor_line_in_fptypes - SOALEN);
					}

					xfOffs_xodd[y2][ind] = X + x + 2;
					xfOffs_xn_xodd[y2][ind] = X + x + 2;
					if(x >= SOALEN - 2) {
						xfOffs_xodd[y2][ind] += (spinor_line_in_fptypes - SOALEN);
						if(local_dir[0]) {
							xfOffs_xn_xodd[y2][ind] += (spinor_line_in_fptypes - SOALEN - nvecs * spinor_line_in_fptypes);
						}
						else {
							xfOffs_xn_xodd[y2][ind] -= SOALEN; // This lane is disabled, just set it within same cache line
							xfmask_xn_xodd[y2] &= ~(1 << ind); // reset the ind bit in the mask
						}
					}

					ybOffs_y0[ind] = X - 3*nvecs*spinor_line_in_fptypes + x; // previous y-neighbor site offsets
					ybOffs_y2[ind] = X - 3*nvecs*spinor_line_in_fptypes + x; // previous y-neighbor site offsets
					if(y == 0) {
						if(local_dir[1]) {
							ybOffs_y0[ind] += Ny*nvecs*spinor_line_in_fptypes;
							ybOffs_y2[ind] += Ny*nvecs*spinor_line_in_fptypes;
						} 
						else {
							ybOffs_y0[ind] = X + x; // This lane is disabled, just set it within same cache line
							ybOffs_y2[ind] = X + x; // This lane is disabled, just set it within same cache line
							ybmask_y0 &= ~(1 << ind); // reset the ind bit in the mask 
							ybmask_y2 &= ~(1 << ind); // reset the ind bit in the mask 
						}
					}
					else if(y <= 2) {
						if(local_dir[1]) {
							ybOffs_y0[ind] += Ny*nvecs*spinor_line_in_fptypes;
						} 
						else {
							ybOffs_y0[ind] = X + x; // This lane is disabled, just set it within same cache line
							ybmask_y0 &= ~(1 << ind); // reset the ind bit in the mask 
						}
					}
					ybOffs_yn0[ind] = X - 3*nvecs*spinor_line_in_fptypes + x; // previous y-neighbor site offsets
					yfOffs_yn[ind] = X +  3*nvecs*spinor_line_in_fptypes + x; // next y-neighbor site offsets
					yfOffs_yn2[ind] = X +  3*nvecs*spinor_line_in_fptypes + x; // next y-neighbor site offsets
					if(y == nyg - 1) {
						if(local_dir[1]) {
							yfOffs_yn[ind] -= Ny*nvecs*spinor_line_in_fptypes;
							yfOffs_yn2[ind] -= Ny*nvecs*spinor_line_in_fptypes;
						} 
						else {
							yfOffs_yn[ind] = X + x; // This lane is disabled, just set it within same cache line
							yfOffs_yn2[ind] = X + x; // This lane is disabled, just set it within same cache line
							yfmask_yn &= ~(1 << ind); // reset the ind bit in the mask 
							yfmask_yn2 &= ~(1 << ind); // reset the ind bit in the mask 
						}
					}
					else if(y >= nyg - 3) {
						if(local_dir[1]) {
							yfOffs_yn[ind] -= Ny*nvecs*spinor_line_in_fptypes;
						} 
						else {
							yfOffs_yn[ind] = X + x; // This lane is disabled, just set it within same cache line
							yfmask_yn &= ~(1 << ind); // reset the ind bit in the mask 
						}
					}
					yfOffs_ynn[ind] = X +  3*nvecs*spinor_line_in_fptypes + x; // next y-neighbor site offsets
					offs[ind] = X + x;  // site offsets for z & t neighbors
					ind++;
				}
			}
			if(nyg == 1) {
				if(!local_dir[1]) ybmask_y0 = yfmask_yn = 0;
			}
		}
		int* 
        getXbOffs(int x, int xodd) 
        {
            return (x == 0 ? xbOffs_x0_xodd[xodd] : xbOffs_xodd[xodd]);
        }
		
        int*
        getXfOffs(int x, int xodd) 
        {
            return (x == nvecs-1 ? xfOffs_xn_xodd[xodd] : xfOffs_xodd[xodd]);
        }
        
		int*
        getYbOffs(int y) 
        {
            int *ret;
			if(nyg == 1 || nyg == 4) 
				ret = (y < 3 ? ybOffs_y0 : ybOffs_yn0);
			else if(nyg == 2) 
				ret = (y >= 4 ? ybOffs_yn0 : (y == 0 ? ybOffs_y0 : ybOffs_y2));
            return ret;
		}
        
		int*
        getYfOffs(int y) 
        {
            int *ret;
			if(nyg == 1 || nyg == 4) 
				ret =  (y >= Ny - 3 ? yfOffs_yn : yfOffs_ynn);
			else if(nyg == 2) 
				ret = (y < Ny - 5 ? yfOffs_ynn : (y == Ny - 2 ? yfOffs_yn : yfOffs_yn2));
            return ret;
		}
        
		int*
        getOffs() 
        {
            return offs;
        }

		unsigned int 
        getXbMask(int x, int xodd) 
        {
            return (x == 0 ? xbmask_x0_xodd[xodd] : -1);
        }
        
		unsigned int 
        getXfMask(int x, int xodd) 
        {
            return (x == nvecs-1 ? xfmask_xn_xodd[xodd] : -1);
        }
        
		unsigned int 
        getYbMask(int y) 
        {
            unsigned int ret;
			if(nyg == 1 || nyg == 4) 
				ret =  (y < 3 ? ybmask_y0 : -1);
			else if(nyg == 2) 
				ret = (y >= 4 ? -1 : (y == 0 ? ybmask_y0 : ybmask_y2));
            return ret;
		}
        
		unsigned int 
        getYfMask(int y) 
        {
            unsigned int ret;
            if(nyg == 1 || nyg == 4) 
				ret =  (y >= Ny - 3 ? yfmask_yn : -1);
			else if(nyg == 2) 
				ret = (y < Ny - 5 ? -1 :(y == Ny - 2 ? yfmask_yn : yfmask_yn2));
            return ret;
		}
};

class GaugeOffsets
{
	private:
		__declspec(align(64)) int gOffs[VECLEN];

	public:
		GaugeOffsets(int nvecs, std::size_t gs) 
        {

			int gauge_line_in_fptypes = gs/sizeof(fptype);
			int nyg = VECLEN/SOALEN;

			for(int y = 0; y < nyg; y++) {
				// Various indexing things
				int ind = y*SOALEN;
				for(int x = 0; x < SOALEN; x++) {
#ifndef USE_PACKED_GAUGES
					gOffs[ind] = nvecs * y * gauge_line_in_fptypes + x; // offsets for gauges
#else
					gOffs[ind] = ind; // this not used really
#endif
					ind++;
				}
			}
		}
        
		int*
        getOffs() 
        {
            return gOffs;
        }
};

#endif // _OFFS_AND_MASKS_H_



#ifndef _UTILS_H_
#define _UTILS_H_

#if (ARCH == sse)

inline __m128d _mm_int2mask_pd(unsigned int msk) {
	static __m128d allOne = _mm_set1_pd(-1.0);
	__m128d ret = _mm_setzero_pd();
	msk = msk & 0x03;

	switch (msk) {
		case 0: ret = _mm_blend_pd(ret, allOne, 0); break;
		case 1: ret = _mm_blend_pd(ret, allOne, 1); break;
		case 2: ret = _mm_blend_pd(ret, allOne, 2); break;
		case 3: ret = _mm_blend_pd(ret, allOne, 3); break;
	}
	return ret;
}

inline __m128 _mm_int2mask_ps(unsigned int msk) {
	static __m128 allOne = _mm_set1_ps(-1.0);
	__m128 ret = _mm_setzero_ps();
	msk = msk & 0x0F;

	switch (msk) {
		case 0: ret = _mm_blend_ps(ret, allOne, 0); break;
		case 1: ret = _mm_blend_ps(ret, allOne, 1); break;
		case 2: ret = _mm_blend_ps(ret, allOne, 2); break;
		case 3: ret = _mm_blend_ps(ret, allOne, 3); break;
		case 4: ret = _mm_blend_ps(ret, allOne, 4); break;
		case 5: ret = _mm_blend_ps(ret, allOne, 5); break;
		case 6: ret = _mm_blend_ps(ret, allOne, 6); break;
		case 7: ret = _mm_blend_ps(ret, allOne, 7); break;
		case 8: ret = _mm_blend_ps(ret, allOne, 8); break;
		case 9: ret = _mm_blend_ps(ret, allOne, 9); break;
		case 10: ret = _mm_blend_ps(ret, allOne, 10); break;
		case 11: ret = _mm_blend_ps(ret, allOne, 11); break;
		case 12: ret = _mm_blend_ps(ret, allOne, 12); break;
		case 13: ret = _mm_blend_ps(ret, allOne, 13); break;
		case 14: ret = _mm_blend_ps(ret, allOne, 14); break;
		case 15: ret = _mm_blend_ps(ret, allOne, 15); break;
	}
	return ret;
}
#endif

#if (ARCH == avx) || (ARCH == avx2)

inline __m256d _mm256_int2mask_pd(unsigned int msk) {
	static __m256d allOne = _mm256_set1_pd(-1.0);
	__m256d ret = _mm256_setzero_pd();
	msk = msk & 0x0F;

	switch (msk) {
		case 0: ret = _mm256_blend_pd(ret, allOne, 0); break;
		case 1: ret = _mm256_blend_pd(ret, allOne, 1); break;
		case 2: ret = _mm256_blend_pd(ret, allOne, 2); break;
		case 3: ret = _mm256_blend_pd(ret, allOne, 3); break;
		case 4: ret = _mm256_blend_pd(ret, allOne, 4); break;
		case 5: ret = _mm256_blend_pd(ret, allOne, 5); break;
		case 6: ret = _mm256_blend_pd(ret, allOne, 6); break;
		case 7: ret = _mm256_blend_pd(ret, allOne, 7); break;
		case 8: ret = _mm256_blend_pd(ret, allOne, 8); break;
		case 9: ret = _mm256_blend_pd(ret, allOne, 9); break;
		case 10: ret = _mm256_blend_pd(ret, allOne, 10); break;
		case 11: ret = _mm256_blend_pd(ret, allOne, 11); break;
		case 12: ret = _mm256_blend_pd(ret, allOne, 12); break;
		case 13: ret = _mm256_blend_pd(ret, allOne, 13); break;
		case 14: ret = _mm256_blend_pd(ret, allOne, 14); break;
		case 15: ret = _mm256_blend_pd(ret, allOne, 15); break;
	}
	return ret;
}


inline __m256 _mm256_int2mask_ps(unsigned int msk) {
	static __m256 allOne = _mm256_set1_ps(-1.0);
	__m256 ret = _mm256_setzero_ps();
	msk = msk & 0x0FF;

	switch (msk) {
                case 0: ret = _mm256_blend_ps(ret, allOne, 0); break;
                case 1: ret = _mm256_blend_ps(ret, allOne, 1); break;
                case 2: ret = _mm256_blend_ps(ret, allOne, 2); break;
                case 3: ret = _mm256_blend_ps(ret, allOne, 3); break;
                case 4: ret = _mm256_blend_ps(ret, allOne, 4); break;
                case 5: ret = _mm256_blend_ps(ret, allOne, 5); break;
                case 6: ret = _mm256_blend_ps(ret, allOne, 6); break;
                case 7: ret = _mm256_blend_ps(ret, allOne, 7); break;
                case 8: ret = _mm256_blend_ps(ret, allOne, 8); break;
                case 9: ret = _mm256_blend_ps(ret, allOne, 9); break;
                case 10: ret = _mm256_blend_ps(ret, allOne, 10); break;
                case 11: ret = _mm256_blend_ps(ret, allOne, 11); break;
                case 12: ret = _mm256_blend_ps(ret, allOne, 12); break;
                case 13: ret = _mm256_blend_ps(ret, allOne, 13); break;
                case 14: ret = _mm256_blend_ps(ret, allOne, 14); break;
                case 15: ret = _mm256_blend_ps(ret, allOne, 15); break;
                case 16: ret = _mm256_blend_ps(ret, allOne, 16); break;
                case 17: ret = _mm256_blend_ps(ret, allOne, 17); break;
                case 18: ret = _mm256_blend_ps(ret, allOne, 18); break;
                case 19: ret = _mm256_blend_ps(ret, allOne, 19); break;
                case 20: ret = _mm256_blend_ps(ret, allOne, 20); break;
                case 21: ret = _mm256_blend_ps(ret, allOne, 21); break;
                case 22: ret = _mm256_blend_ps(ret, allOne, 22); break;
                case 23: ret = _mm256_blend_ps(ret, allOne, 23); break;
                case 24: ret = _mm256_blend_ps(ret, allOne, 24); break;
                case 25: ret = _mm256_blend_ps(ret, allOne, 25); break;
                case 26: ret = _mm256_blend_ps(ret, allOne, 26); break;
                case 27: ret = _mm256_blend_ps(ret, allOne, 27); break;
                case 28: ret = _mm256_blend_ps(ret, allOne, 28); break;
                case 29: ret = _mm256_blend_ps(ret, allOne, 29); break;
                case 30: ret = _mm256_blend_ps(ret, allOne, 30); break;
                case 31: ret = _mm256_blend_ps(ret, allOne, 31); break;
                case 32: ret = _mm256_blend_ps(ret, allOne, 32); break;
                case 33: ret = _mm256_blend_ps(ret, allOne, 33); break;
                case 34: ret = _mm256_blend_ps(ret, allOne, 34); break;
                case 35: ret = _mm256_blend_ps(ret, allOne, 35); break;
                case 36: ret = _mm256_blend_ps(ret, allOne, 36); break;
                case 37: ret = _mm256_blend_ps(ret, allOne, 37); break;
                case 38: ret = _mm256_blend_ps(ret, allOne, 38); break;
                case 39: ret = _mm256_blend_ps(ret, allOne, 39); break;
                case 40: ret = _mm256_blend_ps(ret, allOne, 40); break;
                case 41: ret = _mm256_blend_ps(ret, allOne, 41); break;
                case 42: ret = _mm256_blend_ps(ret, allOne, 42); break;
                case 43: ret = _mm256_blend_ps(ret, allOne, 43); break;
                case 44: ret = _mm256_blend_ps(ret, allOne, 44); break;
                case 45: ret = _mm256_blend_ps(ret, allOne, 45); break;
                case 46: ret = _mm256_blend_ps(ret, allOne, 46); break;
                case 47: ret = _mm256_blend_ps(ret, allOne, 47); break;
                case 48: ret = _mm256_blend_ps(ret, allOne, 48); break;
                case 49: ret = _mm256_blend_ps(ret, allOne, 49); break;
                case 50: ret = _mm256_blend_ps(ret, allOne, 50); break;
                case 51: ret = _mm256_blend_ps(ret, allOne, 51); break;
                case 52: ret = _mm256_blend_ps(ret, allOne, 52); break;
                case 53: ret = _mm256_blend_ps(ret, allOne, 53); break;
                case 54: ret = _mm256_blend_ps(ret, allOne, 54); break;
                case 55: ret = _mm256_blend_ps(ret, allOne, 55); break;
                case 56: ret = _mm256_blend_ps(ret, allOne, 56); break;
                case 57: ret = _mm256_blend_ps(ret, allOne, 57); break;
                case 58: ret = _mm256_blend_ps(ret, allOne, 58); break;
                case 59: ret = _mm256_blend_ps(ret, allOne, 59); break;
                case 60: ret = _mm256_blend_ps(ret, allOne, 60); break;
                case 61: ret = _mm256_blend_ps(ret, allOne, 61); break;
                case 62: ret = _mm256_blend_ps(ret, allOne, 62); break;
                case 63: ret = _mm256_blend_ps(ret, allOne, 63); break;
                case 64: ret = _mm256_blend_ps(ret, allOne, 64); break;
                case 65: ret = _mm256_blend_ps(ret, allOne, 65); break;
                case 66: ret = _mm256_blend_ps(ret, allOne, 66); break;
                case 67: ret = _mm256_blend_ps(ret, allOne, 67); break;
                case 68: ret = _mm256_blend_ps(ret, allOne, 68); break;
                case 69: ret = _mm256_blend_ps(ret, allOne, 69); break;
                case 70: ret = _mm256_blend_ps(ret, allOne, 70); break;
                case 71: ret = _mm256_blend_ps(ret, allOne, 71); break;
                case 72: ret = _mm256_blend_ps(ret, allOne, 72); break;
                case 73: ret = _mm256_blend_ps(ret, allOne, 73); break;
                case 74: ret = _mm256_blend_ps(ret, allOne, 74); break;
                case 75: ret = _mm256_blend_ps(ret, allOne, 75); break;
                case 76: ret = _mm256_blend_ps(ret, allOne, 76); break;
                case 77: ret = _mm256_blend_ps(ret, allOne, 77); break;
                case 78: ret = _mm256_blend_ps(ret, allOne, 78); break;
                case 79: ret = _mm256_blend_ps(ret, allOne, 79); break;
                case 80: ret = _mm256_blend_ps(ret, allOne, 80); break;
                case 81: ret = _mm256_blend_ps(ret, allOne, 81); break;
                case 82: ret = _mm256_blend_ps(ret, allOne, 82); break;
                case 83: ret = _mm256_blend_ps(ret, allOne, 83); break;
                case 84: ret = _mm256_blend_ps(ret, allOne, 84); break;
                case 85: ret = _mm256_blend_ps(ret, allOne, 85); break;
                case 86: ret = _mm256_blend_ps(ret, allOne, 86); break;
                case 87: ret = _mm256_blend_ps(ret, allOne, 87); break;
                case 88: ret = _mm256_blend_ps(ret, allOne, 88); break;
                case 89: ret = _mm256_blend_ps(ret, allOne, 89); break;
                case 90: ret = _mm256_blend_ps(ret, allOne, 90); break;
                case 91: ret = _mm256_blend_ps(ret, allOne, 91); break;
                case 92: ret = _mm256_blend_ps(ret, allOne, 92); break;
                case 93: ret = _mm256_blend_ps(ret, allOne, 93); break;
                case 94: ret = _mm256_blend_ps(ret, allOne, 94); break;
                case 95: ret = _mm256_blend_ps(ret, allOne, 95); break;
                case 96: ret = _mm256_blend_ps(ret, allOne, 96); break;
                case 97: ret = _mm256_blend_ps(ret, allOne, 97); break;
                case 98: ret = _mm256_blend_ps(ret, allOne, 98); break;
                case 99: ret = _mm256_blend_ps(ret, allOne, 99); break;
                case 100: ret = _mm256_blend_ps(ret, allOne, 100); break;
                case 101: ret = _mm256_blend_ps(ret, allOne, 101); break;
                case 102: ret = _mm256_blend_ps(ret, allOne, 102); break;
                case 103: ret = _mm256_blend_ps(ret, allOne, 103); break;
                case 104: ret = _mm256_blend_ps(ret, allOne, 104); break;
                case 105: ret = _mm256_blend_ps(ret, allOne, 105); break;
                case 106: ret = _mm256_blend_ps(ret, allOne, 106); break;
                case 107: ret = _mm256_blend_ps(ret, allOne, 107); break;
                case 108: ret = _mm256_blend_ps(ret, allOne, 108); break;
                case 109: ret = _mm256_blend_ps(ret, allOne, 109); break;
                case 110: ret = _mm256_blend_ps(ret, allOne, 110); break;
                case 111: ret = _mm256_blend_ps(ret, allOne, 111); break;
                case 112: ret = _mm256_blend_ps(ret, allOne, 112); break;
                case 113: ret = _mm256_blend_ps(ret, allOne, 113); break;
                case 114: ret = _mm256_blend_ps(ret, allOne, 114); break;
                case 115: ret = _mm256_blend_ps(ret, allOne, 115); break;
                case 116: ret = _mm256_blend_ps(ret, allOne, 116); break;
                case 117: ret = _mm256_blend_ps(ret, allOne, 117); break;
                case 118: ret = _mm256_blend_ps(ret, allOne, 118); break;
                case 119: ret = _mm256_blend_ps(ret, allOne, 119); break;
                case 120: ret = _mm256_blend_ps(ret, allOne, 120); break;
                case 121: ret = _mm256_blend_ps(ret, allOne, 121); break;
                case 122: ret = _mm256_blend_ps(ret, allOne, 122); break;
                case 123: ret = _mm256_blend_ps(ret, allOne, 123); break;
                case 124: ret = _mm256_blend_ps(ret, allOne, 124); break;
                case 125: ret = _mm256_blend_ps(ret, allOne, 125); break;
                case 126: ret = _mm256_blend_ps(ret, allOne, 126); break;
                case 127: ret = _mm256_blend_ps(ret, allOne, 127); break;
                case 128: ret = _mm256_blend_ps(ret, allOne, 128); break;
                case 129: ret = _mm256_blend_ps(ret, allOne, 129); break;
                case 130: ret = _mm256_blend_ps(ret, allOne, 130); break;
                case 131: ret = _mm256_blend_ps(ret, allOne, 131); break;
                case 132: ret = _mm256_blend_ps(ret, allOne, 132); break;
                case 133: ret = _mm256_blend_ps(ret, allOne, 133); break;
                case 134: ret = _mm256_blend_ps(ret, allOne, 134); break;
                case 135: ret = _mm256_blend_ps(ret, allOne, 135); break;
                case 136: ret = _mm256_blend_ps(ret, allOne, 136); break;
                case 137: ret = _mm256_blend_ps(ret, allOne, 137); break;
                case 138: ret = _mm256_blend_ps(ret, allOne, 138); break;
                case 139: ret = _mm256_blend_ps(ret, allOne, 139); break;
                case 140: ret = _mm256_blend_ps(ret, allOne, 140); break;
                case 141: ret = _mm256_blend_ps(ret, allOne, 141); break;
                case 142: ret = _mm256_blend_ps(ret, allOne, 142); break;
                case 143: ret = _mm256_blend_ps(ret, allOne, 143); break;
                case 144: ret = _mm256_blend_ps(ret, allOne, 144); break;
                case 145: ret = _mm256_blend_ps(ret, allOne, 145); break;
                case 146: ret = _mm256_blend_ps(ret, allOne, 146); break;
                case 147: ret = _mm256_blend_ps(ret, allOne, 147); break;
                case 148: ret = _mm256_blend_ps(ret, allOne, 148); break;
                case 149: ret = _mm256_blend_ps(ret, allOne, 149); break;
                case 150: ret = _mm256_blend_ps(ret, allOne, 150); break;
                case 151: ret = _mm256_blend_ps(ret, allOne, 151); break;
                case 152: ret = _mm256_blend_ps(ret, allOne, 152); break;
                case 153: ret = _mm256_blend_ps(ret, allOne, 153); break;
                case 154: ret = _mm256_blend_ps(ret, allOne, 154); break;
                case 155: ret = _mm256_blend_ps(ret, allOne, 155); break;
                case 156: ret = _mm256_blend_ps(ret, allOne, 156); break;
                case 157: ret = _mm256_blend_ps(ret, allOne, 157); break;
                case 158: ret = _mm256_blend_ps(ret, allOne, 158); break;
                case 159: ret = _mm256_blend_ps(ret, allOne, 159); break;
                case 160: ret = _mm256_blend_ps(ret, allOne, 160); break;
                case 161: ret = _mm256_blend_ps(ret, allOne, 161); break;
                case 162: ret = _mm256_blend_ps(ret, allOne, 162); break;
                case 163: ret = _mm256_blend_ps(ret, allOne, 163); break;
                case 164: ret = _mm256_blend_ps(ret, allOne, 164); break;
                case 165: ret = _mm256_blend_ps(ret, allOne, 165); break;
                case 166: ret = _mm256_blend_ps(ret, allOne, 166); break;
                case 167: ret = _mm256_blend_ps(ret, allOne, 167); break;
                case 168: ret = _mm256_blend_ps(ret, allOne, 168); break;
                case 169: ret = _mm256_blend_ps(ret, allOne, 169); break;
                case 170: ret = _mm256_blend_ps(ret, allOne, 170); break;
                case 171: ret = _mm256_blend_ps(ret, allOne, 171); break;
                case 172: ret = _mm256_blend_ps(ret, allOne, 172); break;
                case 173: ret = _mm256_blend_ps(ret, allOne, 173); break;
                case 174: ret = _mm256_blend_ps(ret, allOne, 174); break;
                case 175: ret = _mm256_blend_ps(ret, allOne, 175); break;
                case 176: ret = _mm256_blend_ps(ret, allOne, 176); break;
                case 177: ret = _mm256_blend_ps(ret, allOne, 177); break;
                case 178: ret = _mm256_blend_ps(ret, allOne, 178); break;
                case 179: ret = _mm256_blend_ps(ret, allOne, 179); break;
                case 180: ret = _mm256_blend_ps(ret, allOne, 180); break;
                case 181: ret = _mm256_blend_ps(ret, allOne, 181); break;
                case 182: ret = _mm256_blend_ps(ret, allOne, 182); break;
                case 183: ret = _mm256_blend_ps(ret, allOne, 183); break;
                case 184: ret = _mm256_blend_ps(ret, allOne, 184); break;
                case 185: ret = _mm256_blend_ps(ret, allOne, 185); break;
                case 186: ret = _mm256_blend_ps(ret, allOne, 186); break;
                case 187: ret = _mm256_blend_ps(ret, allOne, 187); break;
                case 188: ret = _mm256_blend_ps(ret, allOne, 188); break;
                case 189: ret = _mm256_blend_ps(ret, allOne, 189); break;
                case 190: ret = _mm256_blend_ps(ret, allOne, 190); break;
                case 191: ret = _mm256_blend_ps(ret, allOne, 191); break;
                case 192: ret = _mm256_blend_ps(ret, allOne, 192); break;
                case 193: ret = _mm256_blend_ps(ret, allOne, 193); break;
                case 194: ret = _mm256_blend_ps(ret, allOne, 194); break;
                case 195: ret = _mm256_blend_ps(ret, allOne, 195); break;
                case 196: ret = _mm256_blend_ps(ret, allOne, 196); break;
                case 197: ret = _mm256_blend_ps(ret, allOne, 197); break;
                case 198: ret = _mm256_blend_ps(ret, allOne, 198); break;
                case 199: ret = _mm256_blend_ps(ret, allOne, 199); break;
                case 200: ret = _mm256_blend_ps(ret, allOne, 200); break;
                case 201: ret = _mm256_blend_ps(ret, allOne, 201); break;
                case 202: ret = _mm256_blend_ps(ret, allOne, 202); break;
                case 203: ret = _mm256_blend_ps(ret, allOne, 203); break;
                case 204: ret = _mm256_blend_ps(ret, allOne, 204); break;
                case 205: ret = _mm256_blend_ps(ret, allOne, 205); break;
                case 206: ret = _mm256_blend_ps(ret, allOne, 206); break;
                case 207: ret = _mm256_blend_ps(ret, allOne, 207); break;
                case 208: ret = _mm256_blend_ps(ret, allOne, 208); break;
                case 209: ret = _mm256_blend_ps(ret, allOne, 209); break;
                case 210: ret = _mm256_blend_ps(ret, allOne, 210); break;
                case 211: ret = _mm256_blend_ps(ret, allOne, 211); break;
                case 212: ret = _mm256_blend_ps(ret, allOne, 212); break;
                case 213: ret = _mm256_blend_ps(ret, allOne, 213); break;
                case 214: ret = _mm256_blend_ps(ret, allOne, 214); break;
                case 215: ret = _mm256_blend_ps(ret, allOne, 215); break;
                case 216: ret = _mm256_blend_ps(ret, allOne, 216); break;
                case 217: ret = _mm256_blend_ps(ret, allOne, 217); break;
                case 218: ret = _mm256_blend_ps(ret, allOne, 218); break;
                case 219: ret = _mm256_blend_ps(ret, allOne, 219); break;
                case 220: ret = _mm256_blend_ps(ret, allOne, 220); break;
                case 221: ret = _mm256_blend_ps(ret, allOne, 221); break;
                case 222: ret = _mm256_blend_ps(ret, allOne, 222); break;
                case 223: ret = _mm256_blend_ps(ret, allOne, 223); break;
                case 224: ret = _mm256_blend_ps(ret, allOne, 224); break;
                case 225: ret = _mm256_blend_ps(ret, allOne, 225); break;
                case 226: ret = _mm256_blend_ps(ret, allOne, 226); break;
                case 227: ret = _mm256_blend_ps(ret, allOne, 227); break;
                case 228: ret = _mm256_blend_ps(ret, allOne, 228); break;
                case 229: ret = _mm256_blend_ps(ret, allOne, 229); break;
                case 230: ret = _mm256_blend_ps(ret, allOne, 230); break;
                case 231: ret = _mm256_blend_ps(ret, allOne, 231); break;
                case 232: ret = _mm256_blend_ps(ret, allOne, 232); break;
                case 233: ret = _mm256_blend_ps(ret, allOne, 233); break;
                case 234: ret = _mm256_blend_ps(ret, allOne, 234); break;
                case 235: ret = _mm256_blend_ps(ret, allOne, 235); break;
                case 236: ret = _mm256_blend_ps(ret, allOne, 236); break;
                case 237: ret = _mm256_blend_ps(ret, allOne, 237); break;
                case 238: ret = _mm256_blend_ps(ret, allOne, 238); break;
                case 239: ret = _mm256_blend_ps(ret, allOne, 239); break;
                case 240: ret = _mm256_blend_ps(ret, allOne, 240); break;
                case 241: ret = _mm256_blend_ps(ret, allOne, 241); break;
                case 242: ret = _mm256_blend_ps(ret, allOne, 242); break;
                case 243: ret = _mm256_blend_ps(ret, allOne, 243); break;
                case 244: ret = _mm256_blend_ps(ret, allOne, 244); break;
                case 245: ret = _mm256_blend_ps(ret, allOne, 245); break;
                case 246: ret = _mm256_blend_ps(ret, allOne, 246); break;
                case 247: ret = _mm256_blend_ps(ret, allOne, 247); break;
                case 248: ret = _mm256_blend_ps(ret, allOne, 248); break;
                case 249: ret = _mm256_blend_ps(ret, allOne, 249); break;
                case 250: ret = _mm256_blend_ps(ret, allOne, 250); break;
                case 251: ret = _mm256_blend_ps(ret, allOne, 251); break;
                case 252: ret = _mm256_blend_ps(ret, allOne, 252); break;
                case 253: ret = _mm256_blend_ps(ret, allOne, 253); break;
                case 254: ret = _mm256_blend_ps(ret, allOne, 254); break;
                case 255: ret = _mm256_blend_ps(ret, allOne, 255); break;
	}
	return ret;
}

#endif
#endif

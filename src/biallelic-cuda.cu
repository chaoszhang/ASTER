#include "biallelic-cuda.cuh"
#include <cuda_runtime_api.h>
#include <cuda.h>

__constant__ short TAPE[MAX_TAPE_SIZE * 4];

__global__ void updateBatch(GroupIBatch *GI, SiteBatch *S, HotBatch *H, int n, int x, int y, int i){
	int tx = threadIdx.x;
    int bx = blockIdx.x;

	int len = GI[bx][1] - GI[bx][0];
	int start = n * GI[bx][0] + i * len;
	
	for (int j = 0; j < len; j++){
		for (int p = 0; p < 4; p++){
			if (y != -1) S[GI[bx][0] + j][p][y][tx] -= ((H[start + j][p] >> tx) & 1);
			if (x != -1) S[GI[bx][0] + j][p][x][tx] += ((H[start + j][p] >> tx) & 1);
		}
	}
}

void cudaUpdate(GroupIBatch *GI, SiteBatch *S, HotBatch *H, int n, int x, int y, int i, int kc, int b){
	updateBatch<<<kc, b>>>(GI, S, H, n, x, y, i);
}

__global__ void scoreBatch(GroupIBatch *GI, GroupFBatch *GF, SiteBatch *S, double *BS){
	int tx = threadIdx.x;
    int bx = blockIdx.x;

	__shared__ int sXX[7][BATCH], sYZ[7][BATCH], sXpXq[7][BATCH], sYpZq[7][BATCH];
	__shared__ int sYY[7][BATCH], sXZ[7][BATCH], sYpYq[7][BATCH], sXpZq[7][BATCH];
	__shared__ int sZZ[7][BATCH], sXY[7][BATCH], sZpZq[7][BATCH], sXpYq[7][BATCH];
	__shared__ long long m[7][BATCH], mpq[7][BATCH], mp2q2[7][BATCH];
	
	for (int p = 0; p < 7; p++){
		sXX[p][tx] = 0; sYZ[p][tx] = 0; sXpXq[p][tx] = 0; sYpZq[p][tx] = 0;
		sYY[p][tx] = 0; sXZ[p][tx] = 0; sYpYq[p][tx] = 0; sXpZq[p][tx] = 0;
		sZZ[p][tx] = 0; sXY[p][tx] = 0; sZpZq[p][tx] = 0; sXpYq[p][tx] = 0;
		m[p][tx] = 0; mpq[p][tx] = 0; mp2q2[p][tx] = 0;
	}
	for (int j = GI[bx][0], jMax = GI[bx][1]; j < jMax; j++){
		int x = 0, y = 0, z = 0;
		short (&Xh)[4][3][BATCH] = S[j];
		for (int p = 0; p < 4; p++){
			x += Xh[p][0][tx];
			y += Xh[p][1][tx];
			z += Xh[p][2][tx];
		}
		for (int p = 0; p < 7; p++){
			int xp = (p < 4) ? Xh[p][0][tx] : Xh[0][0][tx] + Xh[p-3][0][tx];
			int yp = (p < 4) ? Xh[p][1][tx] : Xh[0][1][tx] + Xh[p-3][1][tx];
			int zp = (p < 4) ? Xh[p][2][tx] : Xh[0][2][tx] + Xh[p-3][2][tx];
			const int xq = x - xp, yq = y - yp, zq = z - zp;
			
			const int xx = x * (x - 1), xpxq = xp * xq, yz = y * z, ypzq = yp * zq + zp * yq;
			const int yy = y * (y - 1), ypyq = yp * yq, xz = x * z, xpzq = xp * zq + zp * xq;
			const int zz = z * (z - 1), zpzq = zp * zq, xy = x * y, xpyq = xp * yq + yp * xq;
			
			sXX[p][tx] += xx; sXpXq[p][tx] += xpxq; sYZ[p][tx] += yz; sYpZq[p][tx] += ypzq;
			sYY[p][tx] += yy; sYpYq[p][tx] += ypyq; sXZ[p][tx] += xz; sXpZq[p][tx] += xpzq;
			sZZ[p][tx] += zz; sZpZq[p][tx] += zpzq; sXY[p][tx] += xy; sXpYq[p][tx] += xpyq;
			
			const long long _xx = xx, _xpxq = xpxq, _yz = yz, _ypzq = ypzq;
			const long long _yy = yy, _ypyq = ypyq, _xz = xz, _xpzq = xpzq;
			const long long _zz = zz, _zpzq = zpzq, _xy = xy, _xpyq = xpyq;

			m[p][tx] += _xpxq * _ypzq + _ypyq * _xpzq + _zpzq * _xpyq;
			mpq[p][tx] -= 2 * (_xpxq * _yz + _ypyq * _xz + _zpzq * _xy)
				 		+ _xx * _ypzq + _yy * _xpzq + _zz * _xpyq;
			mp2q2[p][tx] += 2 * (_xx * _yz + _yy * _xz + _zz * _xy);
		}
	}
	double result = 0.0;
	for (int p = 0; p < 7; p++){
		const double pv = (p < 4) ? GF[bx][p][tx] : GF[bx][0][tx] + GF[bx][p-3][tx]; 
		const double pq = pv * (1 - pv);
		result += (pq * sXX[p][tx] - sXpXq[p][tx]) * (2 * pq * sYZ[p][tx] - sYpZq[p][tx])
				+ (pq * sYY[p][tx] - sYpYq[p][tx]) * (2 * pq * sXZ[p][tx] - sXpZq[p][tx])
				+ (pq * sZZ[p][tx] - sZpZq[p][tx]) * (2 * pq * sXY[p][tx] - sXpYq[p][tx])
				- ((mp2q2[p][tx] * pq + mpq[p][tx]) * pq + m[p][tx]);
	}
	__shared__ double results[BATCH];
	results[tx] = result * GF[bx][4][tx];
	for (int b = 2; b <= BATCH && (tx&(b-1)) == 0; b *= 2){
		results[tx] += results[tx + b/2];
	}
	if (tx == 0) BS[bx] = results[tx];
}

void cudaScore(GroupIBatch *GI, GroupFBatch *GF, SiteBatch* S, double* BS, int kc, int b){
	scoreBatch<<<kc, b>>>(GI, GF, S, BS);
}

void myCudaMalloc(void **p, int size){
	cudaMalloc(p, size);
}

void myCudaMemcpyH2D(void *T, void *S, int size){
	cudaMemcpy(T, S, size, cudaMemcpyHostToDevice);
}

void myCudaMemcpyD2H(void *T, void *S, int size){
	cudaMemcpy(T, S, size, cudaMemcpyDeviceToHost);
}

void myCudaFree(void *p){
	cudaFree(p);
}

void myCudaMemcpyH2C(void *S, int size){
	cudaMemcpyToSymbol(TAPE, S, size);
}

__global__ void cudaWorkBatch(GroupIBatch *GI, GroupFBatch *GF, SiteBatch* S, double* BS, HotBatch *H, int n, int tapeSize){
	int tx = threadIdx.x;
    int bx = blockIdx.x;
	int gs = gridDim.x;

	const int jStart = GI[bx][0], jEnd = GI[bx][1];
	const int len = jEnd - jStart;

	__shared__ int sXX[7][BATCH], sYZ[7][BATCH], sXpXq[7][BATCH], sYpZq[7][BATCH];
	__shared__ int sYY[7][BATCH], sXZ[7][BATCH], sYpYq[7][BATCH], sXpZq[7][BATCH];
	__shared__ int sZZ[7][BATCH], sXY[7][BATCH], sZpZq[7][BATCH], sXpYq[7][BATCH];
	__shared__ long long m[7][BATCH], mpq[7][BATCH], mp2q2[7][BATCH];

	for (int r = 0, s = bx; r < tapeSize; r+=4){
		{
			int x = TAPE[r], y = TAPE[r+1], i = TAPE[r+2];
			int start = n * jStart + i * len;
			for (int j = 0; j < len; j++){
				for (int p = 0; p < 4; p++){
					if (y != -1) S[jStart + j][p][y][tx] -= ((H[start + j][p] >> tx) & 1);
					if (x != -1) S[jStart + j][p][x][tx] += ((H[start + j][p] >> tx) & 1);
				}
			}
		}
		if (TAPE[r+3] > 0){
			for (int p = 0; p < 7; p++){
				sXX[p][tx] = 0; sYZ[p][tx] = 0; sXpXq[p][tx] = 0; sYpZq[p][tx] = 0;
				sYY[p][tx] = 0; sXZ[p][tx] = 0; sYpYq[p][tx] = 0; sXpZq[p][tx] = 0;
				sZZ[p][tx] = 0; sXY[p][tx] = 0; sZpZq[p][tx] = 0; sXpYq[p][tx] = 0;
				m[p][tx] = 0; mpq[p][tx] = 0; mp2q2[p][tx] = 0;
			}
			for (int j = jStart; j < jEnd; j++){
				int x = 0, y = 0, z = 0;
				short (&Xh)[4][3][BATCH] = S[j];
				for (int p = 0; p < 4; p++){
					x += Xh[p][0][tx];
					y += Xh[p][1][tx];
					z += Xh[p][2][tx];
				}
				for (int p = 0; p < 7; p++){
					int xp = (p < 4) ? Xh[p][0][tx] : Xh[0][0][tx] + Xh[p-3][0][tx];
					int yp = (p < 4) ? Xh[p][1][tx] : Xh[0][1][tx] + Xh[p-3][1][tx];
					int zp = (p < 4) ? Xh[p][2][tx] : Xh[0][2][tx] + Xh[p-3][2][tx];
					const int xq = x - xp, yq = y - yp, zq = z - zp;
					
					const int xx = x * (x - 1), xpxq = xp * xq, yz = y * z, ypzq = yp * zq + zp * yq;
					const int yy = y * (y - 1), ypyq = yp * yq, xz = x * z, xpzq = xp * zq + zp * xq;
					const int zz = z * (z - 1), zpzq = zp * zq, xy = x * y, xpyq = xp * yq + yp * xq;
					
					sXX[p][tx] += xx; sXpXq[p][tx] += xpxq; sYZ[p][tx] += yz; sYpZq[p][tx] += ypzq;
					sYY[p][tx] += yy; sYpYq[p][tx] += ypyq; sXZ[p][tx] += xz; sXpZq[p][tx] += xpzq;
					sZZ[p][tx] += zz; sZpZq[p][tx] += zpzq; sXY[p][tx] += xy; sXpYq[p][tx] += xpyq;
					
					const long long _xx = xx, _xpxq = xpxq, _yz = yz, _ypzq = ypzq;
					const long long _yy = yy, _ypyq = ypyq, _xz = xz, _xpzq = xpzq;
					const long long _zz = zz, _zpzq = zpzq, _xy = xy, _xpyq = xpyq;

					m[p][tx] += _xpxq * _ypzq + _ypyq * _xpzq + _zpzq * _xpyq;
					mpq[p][tx] -= 2 * (_xpxq * _yz + _ypyq * _xz + _zpzq * _xy)
								+ _xx * _ypzq + _yy * _xpzq + _zz * _xpyq;
					mp2q2[p][tx] += 2 * (_xx * _yz + _yy * _xz + _zz * _xy);
				}
			}
			double result = 0.0;
			for (int p = 0; p < 7; p++){
				const double pv = (p < 4) ? GF[bx][p][tx] : GF[bx][0][tx] + GF[bx][p-3][tx]; 
				const double pq = pv * (1 - pv);
				result += (pq * sXX[p][tx] - sXpXq[p][tx]) * (2 * pq * sYZ[p][tx] - sYpZq[p][tx])
						+ (pq * sYY[p][tx] - sYpYq[p][tx]) * (2 * pq * sXZ[p][tx] - sXpZq[p][tx])
						+ (pq * sZZ[p][tx] - sZpZq[p][tx]) * (2 * pq * sXY[p][tx] - sXpYq[p][tx])
						- ((mp2q2[p][tx] * pq + mpq[p][tx]) * pq + m[p][tx]);
			}
			__shared__ double results[BATCH];
			results[tx] = result * GF[bx][4][tx];
			for (int b = 2; b <= BATCH && (tx&(b-1)) == 0; b *= 2){
				results[tx] += results[tx + b/2];
			}
			if (tx == 0) BS[s] = results[tx];
			s += gs;
		}
	}
}

void cudaWork(GroupIBatch *GI, GroupFBatch *GF, SiteBatch* S, double* BS, HotBatch *H, int n, int kc, int b, int tapeSize){
	cudaWorkBatch<<<kc, b>>>(GI, GF, S, BS, H, n, tapeSize);
}

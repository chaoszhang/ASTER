#define OBJECTIVE_VERSION "0"

#include<vector>
#include<array>
#include "threadpool.hpp"

using namespace std;

#define CACHE_ALIGNMENT 64

struct TripartitionInitializer{
	vector<int> breakPoints;
	vector<array<vector<bool>, 4> > h;
	vector<array<score_t, 4> > p;
	vector<score_t> w;
	int numThreads = 1;
};

struct ScoreStructure{
	vector<short*> Xh[4][3], X[3], Xp[3];
	vector<score_t*> pq;
};

inline score_t scorePosByLetterColor(const vector<int> &breakPoints, const short *X, const short *Y, const short *Z, 
		const short *Xp, const short *Yp, const short *Zp, const score_t *pq, const vector<score_t> &w, int iMin, int iMax){
	score_t result = 0.0;
	int j = (iMin == 0) ? 0 : breakPoints[iMin - 1];
	for (int i = iMin; i < iMax; i++){
		int sXX = 0, sYZ = 0, sXpXq = 0, sYpZq_sZpYq = 0;
		long long m = 0, mpq = 0, mp2q2 = 0;
		for (int jMax = breakPoints[i]; j < jMax; j++){
			const int x = X[j], y = Y[j], z = Z[j];
			const int xp = Xp[j], yp = Yp[j], zp = Zp[j];
			const int xq = x - xp, yq = y - yp, zq = z - zp;
			sXX += x * (x - 1); sXpXq += xp * xq; 
			sYZ += y * z; sYpZq_sZpYq += yp * zq + zp * yq;
			m += ((long long)(xp * xq)) * ((long long)(yp * zq + zp * yq));
			mpq -= ((long long)(xp * xq)) * ((long long)(2 * y * z)) + ((long long)(x * (x - 1))) * ((long long)(yp * zq + zp * yq));
			mp2q2 += ((long long)(x * (x - 1))) * ((long long)(2 * y * z));
		}
 		result += ((pq[i] * sXX - sXpXq) * (2 * pq[i] * sYZ - sYpZq_sZpYq) - ((mp2q2 * pq[i] + mpq) * pq[i] + m)) * w[i];
	}
	return result; 
}

inline score_t scorePosByLetter(int t, const vector<int> &breakPoints, const vector<short*> (&X)[3], const vector<short*> (&Xp)[3], const score_t *pq, const vector<score_t> &w, int iMin, int iMax){
	score_t result = scorePosByLetterColor(breakPoints, X[0][t], X[1][t], X[2][t], Xp[0][t], Xp[1][t], Xp[2][t], pq, w, iMin, iMax)
		+ scorePosByLetterColor(breakPoints, X[1][t], X[2][t], X[0][t], Xp[1][t], Xp[2][t], Xp[0][t], pq, w, iMin, iMax)
		+ scorePosByLetterColor(breakPoints, X[2][t], X[0][t], X[1][t], Xp[2][t], Xp[0][t], Xp[1][t], pq, w, iMin, iMax);
	return result;
}

inline score_t scorePos(int t, ScoreStructure &SS, const TripartitionInitializer &TI){
	score_t result = 0.0;
	int iMin = t * TI.breakPoints.size() / TI.numThreads;
	int iMax = (t + 1) * TI.breakPoints.size() / TI.numThreads;
	int jMin = (iMin == 0) ? 0 : TI.breakPoints[iMin - 1];
	int jMax = TI.breakPoints[iMax - 1];
	for (int p = 0; p < 4; p++){
		for (int i = iMin; i < iMax; i++){
			SS.pq[t][i] = TI.p[i][p] * (1 - TI.p[i][p]);
		}
		result += scorePosByLetter(t, TI.breakPoints, SS.X, SS.Xh[p], SS.pq[t], TI.w, iMin, iMax);
	}
	for (int p = 1; p < 4; p++){
		for (int i = iMin; i < iMax; i++){
			SS.pq[t][i] = (TI.p[i][p] + TI.p[i][0]) * (1 - TI.p[i][p] - TI.p[i][0]);
		}
		for (int c = 0; c < 3; c++){
			for (int j = jMin; j < jMax; j++){
				SS.Xp[c][t][j] = SS.Xh[p][c][t][j] + SS.Xh[0][c][t][j];
			}
		}
		result += scorePosByLetter(t, TI.breakPoints, SS.X, SS.Xp, SS.pq[t], TI.w, iMin, iMax);
	}
	return result;
}

struct Tripartition{
	const TripartitionInitializer &TI;
	ScoreStructure SS;
	
	bool DeletedFlag = false;
	vector<char*> color;
	
	Tripartition(const TripartitionInitializer &init): TI(init){
		for (int t = 0; t < TI.numThreads; t++){
			color.push_back((char*) aligned_alloc(CACHE_ALIGNMENT, TI.h.size() * sizeof(char)));
			for (int i = 0; i < TI.h.size(); i++) color[t][i] = -1;
			for (int p = 0; p < 4; p++){
				for (int c = 0; c < 3; c++){
					SS.Xh[p][c].push_back((short*) aligned_alloc(CACHE_ALIGNMENT, TI.breakPoints.back() * sizeof(short)));
					for (int i = 0; i < TI.breakPoints.back(); i++) SS.Xh[p][c][t][i] = 0;
				}
			}
			for (int c = 0; c < 3; c++){
				SS.X[c].push_back((short*) aligned_alloc(CACHE_ALIGNMENT, TI.breakPoints.back() * sizeof(short)));
				for (int i = 0; i < TI.breakPoints.back(); i++) SS.X[c][t][i] = 0;
				SS.Xp[c].push_back((short*) aligned_alloc(CACHE_ALIGNMENT, TI.breakPoints.back() * sizeof(short)));
				for (int i = 0; i < TI.breakPoints.back(); i++) SS.Xp[c][t][i] = 0;
			}
			SS.pq.push_back((score_t*) aligned_alloc(CACHE_ALIGNMENT, TI.breakPoints.back() * sizeof(score_t)));
			for (int i = 0; i < TI.breakPoints.back(); i++) SS.pq[t][i] = 0;
		}
	}
	

	~Tripartition(){
		if (DeletedFlag) return;
		for (int t = 0; t < TI.numThreads; t++){
			free(color[t]);
			for (int p = 0; p < 4; p++){
				for (int c = 0; c < 3; c++){
					free(SS.Xh[p][c][t]);
				}
			}
			for (int c = 0; c < 3; c++){
				free(SS.X[c][t]);
				free(SS.Xp[c][t]);
			}
			free(SS.pq[t]);
		}
	}

	Tripartition (Tripartition && T): TI(T.TI), SS(move(T.SS)), color(move(T.color)){
		T.DeletedFlag = true;
	}

    Tripartition& operator= ( Tripartition && ) = delete;
    Tripartition( const Tripartition & ) = delete;
    Tripartition& operator= ( const Tripartition & ) = delete;

	void updatePart(int part, int x, int i){
		int y = color[part][i];
		int iMin = part * TI.breakPoints.size() / TI.numThreads;
		int iMax = (part + 1) * TI.breakPoints.size() / TI.numThreads;
		int jMin = (iMin == 0) ? 0 : TI.breakPoints[iMin - 1];
		int jMax = TI.breakPoints[iMax - 1];
		if (y >= 0){
			for (int p = 0; p < 4; p++){
				for (int j = jMin; j < jMax; j++){
					SS.Xh[p][y][part][j] -= TI.h[i][p][j];
					SS.X[y][part][j] -= TI.h[i][p][j];
				}
			}
		}
		if (x >= 0){
			for (int p = 0; p < 4; p++){
				for (int j = jMin; j < jMax; j++){
					SS.Xh[p][x][part][j] += TI.h[i][p][j];
					SS.X[x][part][j] += TI.h[i][p][j];
				}
			}
		}
		color[part][i] = x;
	}

	void update(int x, int i){
		for (int part = 0; part < TI.numThreads; part++){
			updatePart(part, x, i);
		}
	}
	
	score_t scorePart(int part){
		score_t result = scorePos(part, SS, TI);
		return result;
	}

	score_t score(){
		score_t result = 0;
		for (int part = 0; part < TI.numThreads; part++){
			result += scorePos(part, SS, TI);
		}
		return result;
	}
};

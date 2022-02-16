#include<vector>
#include<array>
#include<thread>

using namespace std;

struct TripartitionInitializer{
	vector<int> breakPoints;
	vector<array<vector<bool>, 4> > h;
	vector<array<score_t, 4> > p;
	vector<score_t> w;
};

struct ScoreStructure{
	vector<short> Xh[4][3], X[3], Xp[3];
	vector<score_t> pq;
};

inline score_t scorePosByLetterColor(const vector<int> &breakPoints, const vector<short> &X, const vector<short> &Y, const vector<short> &Z, 
		const vector<short> &Xp, const vector<short> &Yp, const vector<short> &Zp, const vector<score_t> &pq, const vector<score_t> &w){
	score_t result = 0.0;
	for (int i = 0, j = 0; i < breakPoints.size(); i++){
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

inline score_t scorePosByLetter(const vector<int> &breakPoints, const vector<short> (&X)[3], const vector<short> (&Xp)[3], const vector<score_t> &pq, const vector<score_t> &w){
	score_t result = scorePosByLetterColor(breakPoints, X[0], X[1], X[2], Xp[0], Xp[1], Xp[2], pq, w)
		+ scorePosByLetterColor(breakPoints, X[1], X[2], X[0], Xp[1], Xp[2], Xp[0], pq, w)
		+ scorePosByLetterColor(breakPoints, X[2], X[0], X[1], Xp[2], Xp[0], Xp[1], pq, w);
	return result;
}

inline score_t scorePos(ScoreStructure &SS, const TripartitionInitializer &TI){
	score_t result = 0.0;
	for (int p = 0; p < 4; p++){
		for (int i = 0, iMax = TI.breakPoints.size(); i < iMax; i++){
			SS.pq[i] = TI.p[i][p] * (1 - TI.p[i][p]);
		}
		result += scorePosByLetter(TI.breakPoints, SS.X, SS.Xh[p], SS.pq, TI.w);
	}
	for (int p = 1; p < 4; p++){
		for (int i = 0, iMax = TI.breakPoints.size(); i < iMax; i++){
			SS.pq[i] = (TI.p[i][p] + TI.p[i][0]) * (1 - TI.p[i][p] - TI.p[i][0]);
		}
		for (int c = 0; c < 3; c++){
			for (int j = 0, jMax = TI.breakPoints.back(); j < jMax; j++){
				SS.Xp[c][j] = SS.Xh[p][c][j] + SS.Xh[0][c][j];
			}
		}
		result += scorePosByLetter(TI.breakPoints, SS.X, SS.Xp, SS.pq, TI.w);
	}
	return result;
}
/*
inline score_t scorePos(ScoreStructure &SS, const TripartitionInitializer &TI){
	score_t result = 0.0;
	for (int i = 0, j = 0; i < TI.nSections; i++){
		int sXX = 0, sYY = 0, sZZ = 0, sXpXq = 0, sYpYq = 0, sZpZq = 0;
		int sXpYq = 0, sYpXq = 0, sXpZq = 0, sZpXq = 0, sYpZq = 0, sZpYq = 0;
		for (int jMax = TI.breakPoints[i]; j < jMax; j++){
			const int x = SS.X[j], y = SS.Y[j], z = SS.Z[j];
			const int xp = SS.Xp[j], yp = SS.Yp[j], zp = SS.Zp[j];
			const int xq = x - xp, yq = y - yp, zq = z - zp;
			sXX += x * x; sYY += y * y; sZZ += z * z;
			sXpXq += xp * xq; sYpYq += yp * yq; sZpZq += zp * zq;
			sXpYq += xp * yq; sYpXq += yp * xq; sXpZq += xp * zq;
			sZpXq += zp * xq; sYpZq += yp * zq; sZpYq += zp * yq;
		}
		//const long long L = TI.breakPoints[i] - ((i == 0) ? 0 : TI.breakPoints[i-1]);
		const long long _sXX = sXX, _sYY = sYY, _sZZ = sZZ;
		const long long _sXpXq = sXpXq, _sYpYq = sYpYq, _sZpZq = sZpZq;
		const long long _sXpYq = sXpYq, _sYpXq = sYpXq, _sXpZq = sXpZq;
		const long long _sZpXq = sZpXq, _sYpZq = sYpZq, _sZpYq = sZpYq;
		result += ((_sXpXq * (_sYpZq + _sZpYq) + _sYpYq * (_sXpZq + _sZpXq) + _sZpZq * (_sXpYq + _sYpXq)) * 4.0
				- ((_sXX * (_sYpZq + _sZpYq) + _sYY * (_sXpZq + _sZpXq) + _sZZ * (_sXpYq + _sYpXq)) +
				   (_sXpXq * _sYZ + _sYpYq * _sXZ + _sZpZq * _sXY) * 2) * spq
				+ (_sXX * _sYZ + _sYY * _sXZ + _sZZ * _sXY) * 2 * sp2q2);// / L;
	}
	return result;
}
*/

struct Tripartition{
	const TripartitionInitializer &TI;
	ScoreStructure SS;
	vector<int> color;
	
	Tripartition(const TripartitionInitializer &init): TI(init), color(TI.h.size(), -1){
		for (int p = 0; p < 4; p++){
			for (int c = 0; c < 3; c++){
				SS.Xh[p][c].resize(TI.breakPoints.back(), 0);
			}
		}
		for (int c = 0; c < 3; c++){
			SS.X[c].resize(TI.breakPoints.back(), 0);
			SS.Xp[c].resize(TI.breakPoints.back(), 0);
		}
		SS.pq.resize(TI.breakPoints.size(), 0);
	}
	
	void update(int x, int i){
		int y = color[i];
		if (y >= 0){
			for (int p = 0; p < 4; p++){
				for (int j = 0, jMax = TI.breakPoints.back(); j < jMax; j++){
					SS.Xh[p][y][j] -= TI.h[i][p][j];
					SS.X[y][j] -= TI.h[i][p][j];
				}
			}
		}
		if (x >= 0){
			for (int p = 0; p < 4; p++){
				for (int j = 0, jMax = TI.breakPoints.back(); j < jMax; j++){
					SS.Xh[p][x][j] += TI.h[i][p][j];
					SS.X[x][j] += TI.h[i][p][j];
				}
			}
		}
		color[i] = x;
	}
	
	score_t score(){
		score_t result = scorePos(SS, TI);
		return result;
	}
};

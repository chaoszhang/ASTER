#define OBJECTIVE_VERSION "0"

#include<vector>
#include<array>
#include<thread>
#include "threadpool.hpp"

using namespace std;

inline long long XXYY(long long x0, long long x1, long long x2, long long y0, long long y1, long long y2){
	return x0 * (x0 - 1) * y1 * y2 + x1 * (x1 - 1) * y2 * y0 + x2 * (x2 - 1) * y0 * y1
	     + y0 * (y0 - 1) * x1 * x2 + y1 * (y1 - 1) * x2 * x0 + y2 * (y2 - 1) * x0 * x1;
}

inline score_t scorePos(const array<array<int, 4>, 3> cnt, const array<score_t, 4> pi){
//lst = simplify([sABCD(R, R, Y, Y); sABCD(A, A, Y, Y); sABCD(C, C, R, R); sABCD(A, A, C, C)])
//sol = [Pa^2*Pc^2; -Pc^2*(Pa + Pg)^2; -Pa^2*(Pc + Pt)^2; (Pa + Pg)^2*(Pc + Pt)^2]
//lst2 = simplify([sABCD(Y, Y, R, R); sABCD(Y, Y, A, A); sABCD(R, R, C, C); sABCD(C, C, A, A)])

// (Pa^2+Pg^2)*(Pc^2+Pt^2)*sABCD(R, R, Y, Y)
// -Pr^2*(Pc^2+Pt^2)*sABCD(A, A, Y, Y) -Pr^2*(Pc^2+Pt^2)*sABCD(G, G, Y, Y) -(Pa^2+Pg^2)*Py^2*sABCD(C, C, R, R) -(Pa^2+Pg^2)*Py^2*sABCD(T, T, R, R)
// +Pr^2*Py^2*sABCD(A, A, C, C) +Pr^2*Py^2*sABCD(G, G, C, C) +Pr^2*Py^2*sABCD(A, A, T, T) +Pr^2*Py^2*sABCD(G, G, T, T)

	const score_t A = pi[0], C = pi[1], G = pi[2], T = pi[3];
	const score_t R = A + G, Y = C + T, R2 = A * A + G * G, Y2 = C * C + T * T;
	const long long a0 = cnt[0][0], c0 = cnt[0][1], g0 = cnt[0][2], t0 = cnt[0][3], r0 = a0 + g0, y0 = c0 + t0;
	const long long a1 = cnt[1][0], c1 = cnt[1][1], g1 = cnt[1][2], t1 = cnt[1][3], r1 = a1 + g1, y1 = c1 + t1;
	const long long a2 = cnt[2][0], c2 = cnt[2][1], g2 = cnt[2][2], t2 = cnt[2][3], r2 = a2 + g2, y2 = c2 + t2;
	
	const long long rryy = XXYY(r0, r1, r2, y0, y1, y2);

	const long long aayy = XXYY(a0, a1, a2, y0, y1, y2);
	const long long ggyy = XXYY(g0, g1, g2, y0, y1, y2);
	const long long rrcc = XXYY(r0, r1, r2, c0, c1, c2);
	const long long rrtt = XXYY(r0, r1, r2, t0, t1, t2);
	
	const long long aacc = XXYY(a0, a1, a2, c0, c1, c2);
	const long long aatt = XXYY(a0, a1, a2, t0, t1, t2);
	const long long ggcc = XXYY(g0, g1, g2, c0, c1, c2);
	const long long ggtt = XXYY(g0, g1, g2, t0, t1, t2);
	
	return rryy * R2 * Y2 - (aayy + ggyy) * (R * R) * Y2 - (rrcc + rrtt) * R2 * (Y * Y)
	     + (aacc + aatt + ggcc + ggtt) * (R * R) * (Y * Y);
}

struct TripartitionInitializer{
	int npos = 0, nThreads = 1;
	vector<score_t> w;
	vector<array<score_t, 4> > pi;
	vector<vector<char> > seq;
	vector<vector<int> > speciesMapper;
};

struct Tripartition{
	const TripartitionInitializer& TI;
	vector<vector<int> > color;
	vector<vector<array<array<int, 4>, 3> > > cnt; // tree * pos * part * letter
	vector<vector<bool> > valid; // thread * tree
	vector<vector<score_t> > scores; // thread * tree
	// i, j, k -> taxon; p, q, r -> pos; x, y, z -> part; a, b, c -> letter
	
	Tripartition(const TripartitionInitializer &init): TI(init), color(init.nThreads, vector<int>(init.seq.size(), -1)), 
		cnt(init.w.size(), vector<array<array<int, 4>, 3> >(init.npos)),
		valid(init.nThreads, vector<bool>(init.w.size())), scores(init.nThreads, vector<score_t>(init.w.size())){}
	
	void updatePart(int part, int x, int i){
		int y = color[part][i];
		int start = TI.npos * part / TI.nThreads, end = TI.npos * (1 + part) / TI.nThreads;
		color[part][i] = x;
		for (int t: TI.speciesMapper[i]){
			valid[part][t] = false;
			for (int p = start; p < end; p++){
				int a = TI.seq[i][p];
				if (y != -1 && a != -1) cnt[t][p][y][a]--;
				if (x != -1 && a != -1) cnt[t][p][x][a]++;
			}
		}
	}
	
	score_t scorePart(int part){
		int start = TI.npos * part / TI.nThreads, end = TI.npos * (1 + part) / TI.nThreads;
		score_t result = 0;
		for (int t = 0; t < TI.w.size(); t++){	
			if (valid[part][t]){
				result += scores[part][t];
				continue;
			}
			score_t temp = 0;
			for (int p = start; p < end; p++){
				temp += scorePos(cnt[t][p], TI.pi[p]);
			}
			scores[part][t] = temp * TI.w[t];
			valid[part][t] = true;
			result += scores[part][t];
		}
		return result;
	}
};

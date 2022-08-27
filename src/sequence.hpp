#define OBJECTIVE_VERSION "0"

#include<vector>
#include<array>
#include<thread>
#include<future>
#include "threadpool.hpp"

//#define G_SUPPORT

using namespace std;

inline long long XXYY(long long x0, long long x1, long long x2, long long y0, long long y1, long long y2){
	return x0 * (x0 - 1) * y1 * y2 + x1 * (x1 - 1) * y2 * y0 + x2 * (x2 - 1) * y0 * y1
	     + y0 * (y0 - 1) * x1 * x2 + y1 * (y1 - 1) * x2 * x0 + y2 * (y2 - 1) * x0 * x1;
}

inline score_t scorePos(const array<array<unsigned short, 4>, 3> cnt, const array<float, 4> pi){
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
	struct Tree{
		int npos = 0, ntree = 0;
		vector<array<float, 4> > pi; // pos * letter
		vector<vector<char> > seq; // taxon * letter
		vector<vector<unsigned short> > speciesMapper; // taxon * tree

		Tree(){}
		Tree(int n): seq(n), speciesMapper(n){}
	};

	int nThreads = 1;
	vector<vector<Tree> > trees;
};

struct Tripartition{
	struct Tree{
		vector<vector<array<array<unsigned short, 4>, 3> > > cnt; // tree * pos * part * letter
		vector<score_t> scores; // tree
		vector<bool> valid; // tree
		
		Tree(){}
		Tree(int ntree, int npos): scores(ntree), valid(ntree),
			cnt(ntree, vector<array<array<unsigned short, 4>, 3> >(npos)){}
	};

	const TripartitionInitializer& TI;
	vector<vector<char> > color;
	vector<vector<Tree> > trees; // gene * bin
	
	// i, j, k -> taxon; p, q, r -> pos; x, y, z -> part; a, b, c -> letter
	
	Tripartition(const TripartitionInitializer &init): TI(init), trees(init.trees.size()),
			color(init.nThreads, vector<char>(init.trees[0][0].seq.size(), -1)){
		for (int i = 0; i < TI.trees.size(); i++){
			for (int j = 0; j < TI.trees[i].size(); j++){
				trees[i].emplace_back(TI.trees[i][j].ntree, TI.trees[i][j].npos);
			}
		}
	}

	void updatePart(int part, int x, int i){
		int start = trees.size() * part / TI.nThreads, end = trees.size() * (1 + part) / TI.nThreads;
		int y = color[part][i];
		color[part][i] = x;
		for (int a = start; a < end; a++){
			for (int b = 0; b < trees[a].size(); b++){
				for (int t: TI.trees[a][b].speciesMapper[i]){
					trees[a][b].valid[t] = false;
					for (int p = 0; p < trees[a][b].cnt[t].size(); p++){
						int c = TI.trees[a][b].seq[i][p];
						if (c == -1) continue;
						if (y != -1) trees[a][b].cnt[t][p][y][c]--;
						if (x != -1) trees[a][b].cnt[t][p][x][c]++;
					}
				}
			}
		}
	}
	
	score_t scorePart(int part){
		int start = trees.size() * part / TI.nThreads, end = trees.size() * (1 + part) / TI.nThreads;
		score_t result = 0;
		for (int a = start; a < end; a++){
			for (int b = 0; b < trees[a].size(); b++){
				for (int t = 0; t < trees[a][b].scores.size(); t++){
					if (trees[a][b].valid[t]){
						result += trees[a][b].scores[t];
						continue;
					}
					score_t temp = 0;
					for (int p = 0; p < trees[a][b].cnt[t].size(); p++){
						temp += scorePos(trees[a][b].cnt[t][p], TI.trees[a][b].pi[p]);
					}
					trees[a][b].scores[t] = temp;
					trees[a][b].valid[t] = true;
					result += trees[a][b].scores[t];
				}
			}
		}
		return result;
	}
};
/*
inline long long supportXXYY(long long x0, long long x1, long long x2, long long x3,
		long long y0, long long y1, long long y2, long long y3){
	return x0 * x1 * y2 * y3 + y0 * y1 * x2 * x3;
}

inline score_t supportPos(const array<array<unsigned short, 4>, 4> cnt, const array<score_t, 4> pi){
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
	const long long a3 = cnt[3][0], c3 = cnt[3][1], g3 = cnt[3][2], t3 = cnt[3][3], r3 = a3 + g3, y3 = c3 + t3;
	
	const long long rryy = supportXXYY(r0, r1, r2, r3, y0, y1, y2, y3);

	const long long aayy = supportXXYY(a0, a1, a2, a3, y0, y1, y2, y3);
	const long long ggyy = supportXXYY(g0, g1, g2, g3, y0, y1, y2, y3);
	const long long rrcc = supportXXYY(r0, r1, r2, r3, c0, c1, c2, c3);
	const long long rrtt = supportXXYY(r0, r1, r2, r3, t0, t1, t2, t3);
	
	const long long aacc = supportXXYY(a0, a1, a2, a3, c0, c1, c2, c3);
	const long long aatt = supportXXYY(a0, a1, a2, a3, t0, t1, t2, t3);
	const long long ggcc = supportXXYY(g0, g1, g2, g3, c0, c1, c2, c3);
	const long long ggtt = supportXXYY(g0, g1, g2, g3, t0, t1, t2, t3);
	
	return rryy * R2 * Y2 - (aayy + ggyy) * (R * R) * Y2 - (rrcc + rrtt) * R2 * (Y * Y)
	     + (aacc + aatt + ggcc + ggtt) * (R * R) * (Y * Y);
}

struct Quadrupartition{
	const TripartitionInitializer& TI;
	vector<vector<char> > color;
	vector<vector<array<array<unsigned short, 4>, 4> > > cnt; // tree * pos * part * letter
	vector<vector<bool> > valid; // thread * tree
	vector<vector<score_t> > scores0, scores1, scores2; // thread * tree
	// i, j, k -> taxon; p, q, r -> pos; x, y, z -> part; a, b, c -> letter
	
	Quadrupartition(const TripartitionInitializer &init): TI(init), color(init.nThreads, vector<char>(init.seq.size(), -1)), 
		cnt(init.w.size(), vector<array<array<unsigned short, 4>, 4> >(init.npos)),
		valid(init.nThreads, vector<bool>(init.w.size())), scores0(init.nThreads, vector<score_t>(init.w.size())),
		scores1(init.nThreads, vector<score_t>(init.w.size())), scores2(init.nThreads, vector<score_t>(init.w.size())){}

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
	
	array<double, 3> scorePart(int part){
		int start = TI.npos * part / TI.nThreads, end = TI.npos * (1 + part) / TI.nThreads;
		array<double, 3> result;
		result[0] = 0;
		result[1] = 0;
		result[2] = 0;
		for (int t = 0; t < TI.w.size(); t++){	
			if (valid[part][t]){
				result[0] += scores0[part][t];
				result[1] += scores1[part][t];
				result[2] += scores2[part][t];
				continue;
			}
			score_t temp0 = 0, temp1 = 0, temp2 = 0;
			for (int p = start; p < end; p++){
				const array<array<unsigned short, 4>, 4> cnt0 = cnt[t][p];
				const array<array<unsigned short, 4>, 4> cnt1 = {cnt0[0], cnt0[2], cnt0[1], cnt0[3]};
				const array<array<unsigned short, 4>, 4> cnt2 = {cnt0[0], cnt0[3], cnt0[1], cnt0[2]};
				temp0 += supportPos(cnt0, TI.pi[p]);
				temp1 += supportPos(cnt1, TI.pi[p]);
				temp2 += supportPos(cnt2, TI.pi[p]);
			}
			scores0[part][t] = temp0 * TI.w[t];
			scores1[part][t] = temp1 * TI.w[t];
			scores2[part][t] = temp2 * TI.w[t];
			valid[part][t] = true;
			result[0] += scores0[part][t];
			result[1] += scores1[part][t];
			result[2] += scores2[part][t];
		}
		return result;
	}

	void update(int x, int i){
		vector<thread> thrds;
		for (int p = 1; p < TI.nThreads; p++) thrds.emplace_back(&Quadrupartition::updatePart, this, p, x, i);
		updatePart(0, x, i);
		for (thread &t: thrds) t.join();
	}
	
	array<double, 3> score(){
		vector<thread> thrds;
		vector<future<array<double, 3> > > temp(TI.nThreads);
		for (int p = 1; p < TI.nThreads; p++){
			temp[p] = async(launch::async, &Quadrupartition::scorePart, this, p);
		}
		array<double, 3> res = scorePart(0);
		for (int p = 1; p < TI.nThreads; p++){
			array<double, 3> t = temp[p].get();
			res[0] += t[0];
			res[1] += t[1];
			res[2] += t[2];
		}
		return res;
	}
};
*/
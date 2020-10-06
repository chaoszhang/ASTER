#include<vector>
#include<array>
#include<thread>

using namespace std;

inline score_t scorePos(const array<array<int, 4>, 3> cnt, const array<score_t, 4> pi){
	const score_t A = pi[0], C = pi[1], G = pi[2], T = pi[3];
	const long long a0 = cnt[0][0], c0 = cnt[0][1], g0 = cnt[0][2], t0 = cnt[0][3];
	const long long a1 = cnt[1][0], c1 = cnt[1][1], g1 = cnt[1][2], t1 = cnt[1][3];
	const long long a2 = cnt[2][0], c2 = cnt[2][1], g2 = cnt[2][2], t2 = cnt[2][3];
	
	const score_t aacc = a0 * (a0 - 1) * c1 * c2 + a1 * (a1 - 1) * c2 * c0 + a2 * (a2 - 1) * c0 * c1;
	const score_t aatt = a0 * (a0 - 1) * t1 * t2 + a1 * (a1 - 1) * t2 * t0 + a2 * (a2 - 1) * t0 * t1;
	const score_t ggcc = g0 * (g0 - 1) * c1 * c2 + g1 * (g1 - 1) * c2 * c0 + g2 * (g2 - 1) * c0 * c1;
	const score_t ggtt = g0 * (g0 - 1) * t1 * t2 + g1 * (g1 - 1) * t2 * t0 + g2 * (g2 - 1) * t0 * t1;
	
	return (aacc * G*G * T*T - aatt * G*G * C*C - ggcc * A*A * T*T + ggtt * A*A * C*C) / (A - G) / (C - T);
}

struct TripartitionInitializer{
	int npos, nThreads = 1;
	vector<array<score_t, 4> > pi;
	vector<vector<int> > seq;
	vector<bool> weight;
};

struct Tripartition{
	int nThreads;
	const int npos;
	vector<array<array<int, 4>, 3> > cnt; // pos * part * letter
	vector<array<int, 4> > total; // pos * letter
	const vector<array<score_t, 4> > pi; // pos * letter
	const vector<vector<int> > seq; // taxon * pos -> letter
	const vector<bool> weight;
	// i, j, k -> taxon; p, q, r -> pos; x, y, z -> part; a, b, c -> letter
	
	Tripartition(const TripartitionInitializer &init): nThreads(init.nThreads), npos(init.npos), seq(init.seq), pi(init.pi), weight(init.weight), cnt(init.npos), total(init.npos){}
	
	void resetWork(int start, int end){
		for (int p = start; p < end; p++){
			for (int a = 0; a < 4; a++){
				cnt[p][0][a] = total[p][a];
				cnt[p][1][a] = 0;
				cnt[p][2][a] = 0;
			}
		}
	}
	
	void addTotalWork(int i, int start, int end){
		for (int p = start; p < end; p++){
			int a = seq[i][p];
			if (a != -1) total[p][a]++;
		}
	}
	
	void rmvTotalWork(int i, int start, int end){
		for (int p = start; p < end; p++){
			int a = seq[i][p];
			if (a != -1) total[p][a]--;
		}
	}
	
	void addWork(int x, int i, int start, int end){
		for (int p = start; p < end; p++){
			int a = seq[i][p];
			if (a != -1) cnt[p][x][a]++;
		}
	}
	
	void rmvWork(int x, int i, int start, int end){
		for (int p = start; p < end; p++){
			int a = seq[i][p];
			if (a != -1) cnt[p][x][a]--;
		}
	}
	
	void scoreWork(score_t &res, int start, int end){
		score_t temp = 0;
		for (int p = start; p < end; p++){
			if (weight[p]) temp += scorePos(cnt[p], pi[p]);
		}
		res = temp;
	}
	
	void reset(){
		vector<thread> thrds;
		for (int t = 1; t < nThreads; t++) thrds.emplace_back(&Tripartition::resetWork, this, npos * t / nThreads, npos * (t + 1) / nThreads);
		resetWork(0, npos / nThreads);
		for (thread &t: thrds) t.join();
	}
	
	void addTotal(int i){
		vector<thread> thrds;
		for (int t = 1; t < nThreads; t++) thrds.emplace_back(&Tripartition::addTotalWork, this, i, npos * t / nThreads, npos * (t + 1) / nThreads);
		addTotalWork(i, 0, npos / nThreads);
		for (thread &t: thrds) t.join();
	}
	
	void rmvTotal(int i){
		vector<thread> thrds;
		for (int t = 1; t < nThreads; t++) thrds.emplace_back(&Tripartition::rmvTotalWork, this, i, npos * t / nThreads, npos * (t + 1) / nThreads);
		rmvTotalWork(i, 0, npos / nThreads);
		for (thread &t: thrds) t.join();
	}
	
	void add(int x, int i){
		vector<thread> thrds;
		for (int t = 1; t < nThreads; t++) thrds.emplace_back(&Tripartition::addWork, this, x, i, npos * t / nThreads, npos * (t + 1) / nThreads);
		addWork(x, i, 0, npos / nThreads);
		for (thread &t: thrds) t.join();
	}
	
	void rmv(int x, int i){
		vector<thread> thrds;
		for (int t = 1; t < nThreads; t++) thrds.emplace_back(&Tripartition::rmvWork, this, x, i, npos * t / nThreads, npos * (t + 1) / nThreads);
		rmvWork(x, i, 0, npos / nThreads);
		for (thread &t: thrds) t.join();
	}
	
	score_t score(){
		vector<score_t> res(nThreads);
		vector<thread> thrds;
		for (int t = 1; t < nThreads; t++) thrds.emplace_back(&Tripartition::scoreWork, this, ref(res[t]), npos * t / nThreads, npos * (t + 1) / nThreads);
		scoreWork(res[0], 0, npos / nThreads);
		for (thread &t: thrds) t.join();
		for (int t = 1; t < nThreads; t++) res[0] += res[t];
		return res[0];
	}
};

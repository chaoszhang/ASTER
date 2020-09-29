#include<vector>
#include<array>
#include<thread>

using namespace std;

inline score_t scorePos(const array<array<int, 4>, 3> cnt, const array<score_t, 4> pi){
	count_t sum[3] = {0, 0, 0};
	for (int x = 0; x < 3; x++){
		for (int a = 0; a < 4; a++){
			sum[x] += cnt[x][a];
		}
	}
	/*
	score_t res = 0;
	for (int a = 0; a < 4; a++){
		for (int b = 0; b < 4; b++){
			if (b == a) continue;
			res += (1 - pi[b]) * ( ((score_t) cnt[0][a]) * (cnt[0][a] - 1) * cnt[1][b] * cnt[2][b] + ((score_t) cnt[0][b]) * (cnt[0][b] - 1) * cnt[1][a] * cnt[2][a]
			                     + ((score_t) cnt[1][a]) * (cnt[1][a] - 1) * cnt[2][b] * cnt[0][b] + ((score_t) cnt[1][b]) * (cnt[1][b] - 1) * cnt[2][a] * cnt[0][a]
			                     + ((score_t) cnt[2][a]) * (cnt[2][a] - 1) * cnt[0][b] * cnt[1][b] + ((score_t) cnt[2][b]) * (cnt[2][b] - 1) * cnt[0][a] * cnt[1][a] )
			       - pi[a] * (pi[a] + pi[b]) * ( ((score_t) cnt[0][b]) * (cnt[0][b] - 1) * sum[1] * sum[2] + ((score_t) sum[0]) * (sum[0] - 1) * cnt[1][b] * cnt[2][b]
                                               + ((score_t) cnt[1][b]) * (cnt[1][b] - 1) * sum[2] * sum[0] + ((score_t) sum[1]) * (sum[1] - 1) * cnt[2][b] * cnt[0][b]
                                               + ((score_t) cnt[2][b]) * (cnt[2][b] - 1) * sum[0] * sum[1] + ((score_t) sum[2]) * (sum[2] - 1) * cnt[0][b] * cnt[1][b] );
		}
	}
	return res;
	*/
	score_t spi = 0;
	score_t p0v0 = 0, p0v1 = 0, p0v2 = 0, p1v0 = 0, p1v1 = 0, p1v2 = 0, p2v0 = 0, p2v1 = 0, p2v2 = 0;
	score_t q0v0 = 0, q0v1 = 0, q0v2 = 0, q1v0 = 0, q1v1 = 0, q1v2 = 0, q2v0 = 0, q2v1 = 0, q2v2 = 0;
	for (int a = 0; a < 4; a++){
		spi += pi[a] * pi[a];
	}
	for (int a = 0; a < 4; a++){
		count_t cnt0 = cnt[0][a], cnt1 = cnt[1][a], cnt2 = cnt[2][a];
		score_t p0 = cnt0 * (cnt0 - 1LL), p1 = cnt1 * (cnt1 - 1LL), p2 = cnt2 * (cnt2 - 1LL);
		score_t q0 = cnt1 * cnt2, q1 = cnt2 * cnt0, q2 = cnt0 * cnt1;
		score_t npi = 1 - pi[a], mpi = spi + pi[a] - 2 * pi[a] * pi[a];
		p0v0 += p0; p0v1 += npi * p0; p0v2 += mpi * p0; q0v0 += q0; q0v1 += npi * q0; q0v2 += mpi * q0;
		p1v0 += p1; p1v1 += npi * p1; p1v2 += mpi * p1; q1v0 += q1; q1v1 += npi * q1; q1v2 += mpi * q1;
		p2v0 += p2; p2v1 += npi * p2; p2v2 += mpi * p2; q2v0 += q2; q2v1 += npi * q2; q2v2 += mpi * q2;
	}
	return p0v0 * q0v1 + p0v1 * q0v0 + p1v0 * q1v1 + p1v1 * q1v0 + p2v0 * q2v1 + p2v1 * q2v0
		- (sum[1] * sum[2] * p0v2 + sum[0] * (sum[0] - 1) * q0v2 + sum[2] * sum[0] * p1v2 + sum[1] * (sum[1] - 1) * q1v2 + sum[0] * sum[1] * p2v2 + sum[2] * (sum[2] - 1) * q2v2);
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
		for (int p = 0; p < npos; p++){
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

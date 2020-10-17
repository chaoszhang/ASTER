#include<vector>
#include<array>
#include<thread>

using namespace std;

inline score_t scorePos(const array<array<int, 4>, 3> cnt, const score_t R){
	const score_t M = 1 - 2 * R + 2 * R * R;
	const long long rr0 = cnt[0][0], ry0 = cnt[0][1], yr0 = cnt[0][2], yy0 = cnt[0][3];
	const long long rr1 = cnt[1][0], ry1 = cnt[1][1], yr1 = cnt[1][2], yy1 = cnt[1][3];
	const long long rr2 = cnt[2][0], ry2 = cnt[2][1], yr2 = cnt[2][2], yy2 = cnt[2][3];
	const long long rn0 = rr0 + ry0, nr0 = rr0 + yr0, yn0 = yr0 + yy0, ny0 = ry0 + yy0;
	const long long rn1 = rr1 + ry1, nr1 = rr1 + yr1, yn1 = yr1 + yy1, ny1 = ry1 + yy1;
	const long long rn2 = rr2 + ry2, nr2 = rr2 + yr2, yn2 = yr2 + yy2, ny2 = ry2 + yy2;
	
	return rn0 * (rn0 - 1) * yn1 * yn2 + nr0 * (nr0 - 1) * ny1 * ny2 - M * (rn0 * (rn0 - 1) * ny1 * ny2 + nr0 * (nr0 - 1) * yn1 * yn2)
	     + rn1 * (rn1 - 1) * yn2 * yn0 + nr1 * (nr1 - 1) * ny2 * ny0 - M * (rn1 * (rn1 - 1) * ny2 * ny0 + nr1 * (nr1 - 1) * yn2 * yn0)
	     + rn2 * (rn2 - 1) * yn0 * yn1 + nr2 * (nr2 - 1) * ny0 * ny1 - M * (rn2 * (rn2 - 1) * ny0 * ny1 + nr2 * (nr2 - 1) * yn0 * yn1);
}

struct TripartitionInitializer{
	int npos, nThreads = 1;
	vector<score_t> pi;
	vector<vector<int> > seq;
	vector<bool> weight;
	vector<score_t> weightHelper;
};

struct Tripartition{
	int nThreads, n = 0;
	const int npos;
	vector<array<array<int, 4>, 3> > cnt; // pos * part * letter
	vector<array<int, 4> > total; // pos * letter
	const vector<score_t> pi; // pos * letter
	const vector<vector<int> > seq; // taxon * pos -> letter
	const vector<bool> weight;
	const vector<score_t> weightHelper;
	// i, j, k -> taxon; p, q, r -> pos; x, y, z -> part; a, b, c -> letter
	
	Tripartition(const TripartitionInitializer &init): nThreads(init.nThreads), npos(init.npos), seq(init.seq), pi(init.pi), weight(init.weight), cnt(init.npos), total(init.npos), weightHelper(init.weightHelper){}
	
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
		score_t temp = 0, tempH = 0;
		for (int p = start; p < end; p++){
			if (weight[p]) {
				temp += scorePos(cnt[p], pi[p]);
				//tempH += scorePosHelper(cnt[p], pi[p]);
			}
		}
		res = temp; // + tempH * weightHelper[n];
	}
	
	void reset(){
		vector<thread> thrds;
		for (int t = 1; t < nThreads; t++) thrds.emplace_back(&Tripartition::resetWork, this, npos * t / nThreads, npos * (t + 1) / nThreads);
		resetWork(0, npos / nThreads);
		for (thread &t: thrds) t.join();
	}
	
	void addTotal(int i){
		n++;
		vector<thread> thrds;
		for (int t = 1; t < nThreads; t++) thrds.emplace_back(&Tripartition::addTotalWork, this, i, npos * t / nThreads, npos * (t + 1) / nThreads);
		addTotalWork(i, 0, npos / nThreads);
		for (thread &t: thrds) t.join();
	}
	
	void rmvTotal(int i){
		n--;
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

#include<vector>
#include<array>
#include<thread>

using namespace std;

inline score_t scorePos(const array<array<int, 4>, 3> cnt, const score_t pi){
	const score_t R = pi, Y = 1 - R;
	const long long r0 = cnt[0][0], p0 = cnt[0][1], q0 = cnt[0][2], y0 = cnt[0][3], n0 = r0 + p0 + q0 + y0;
	const long long r1 = cnt[1][0], p1 = cnt[1][1], q1 = cnt[1][2], y1 = cnt[1][3], n1 = r1 + p1 + q1 + y1;
	const long long r2 = cnt[2][0], p2 = cnt[2][1], q2 = cnt[2][2], y2 = cnt[2][3], n2 = r2 + p2 + q2 + y2;
	
	const score_t rrnn = r0 * (r0 - 1) * n1 * n2 + r1 * (r1 - 1) * n2 * n0 + r2 * (r2 - 1) * n0 * n1;
	const score_t yynn = y0 * (y0 - 1) * n1 * n2 + y1 * (y1 - 1) * n2 * n0 + y2 * (y2 - 1) * n0 * n1;
	const score_t rryy = r0 * (r0 - 1) * y1 * y2 + r1 * (r1 - 1) * y2 * y0 + r2 * (r2 - 1) * y0 * y1;
	const score_t ppqq = p0 * (p0 - 1) * q1 * q2 + p1 * (p1 - 1) * q2 * q0 + p2 * (p2 - 1) * q0 * q1;
	
	const score_t rrpq = p0 * q0 * r1 * r2 + p1 * q1 * r2 * r0 + p2 * q2 * r0 * r1;
	const score_t yypq = p0 * q0 * y1 * y2 + p1 * q1 * y2 * y0 + p2 * q2 * y0 * y1;
	
	return ((R > Y) ? 1 : -1) * (R*R * yynn - Y*Y * rrnn + (R - Y)*(2*rryy - ppqq) + (Y*Y - R*R)*(rryy + 2*rrpq + 2*yypq));
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

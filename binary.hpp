#include<vector>
#include<array>
#include<thread>

using namespace std;

score_t scoreGene(const vector<array<array<int, 2>, 3> > &cnt, const score_t R){
	long long rryy0 = 0, rrnn0 = 0, nnyy0 = 0, nnnn0 = 0;
	long long rryy1 = 0, rrnn1 = 0, nnyy1 = 0, nnnn1 = 0;
	long long rryy2 = 0, rrnn2 = 0, nnyy2 = 0, nnnn2 = 0;
	for (const array<array<int, 2>, 3> &c: cnt){
		const long long r0 = c[0][0], y0 = c[0][1], n0 = r0 + y0;
		const long long r1 = c[1][0], y1 = c[1][1], n1 = r1 + y1;
		const long long r2 = c[2][0], y2 = c[2][1], n2 = r2 + y2;
		rryy0 += r0 * (r0 - 1) * y1 * y2; rrnn0 += r0 * (r0 - 1) * n1 * n2; nnyy0 += n0 * (n0 - 1) * y1 * y2; nnnn0 += n0 * (n0 - 1) * n1 * n2;
		rryy1 += r1 * (r1 - 1) * y2 * y0; rrnn1 += r1 * (r1 - 1) * n2 * n0; nnyy1 += n1 * (n1 - 1) * y2 * y0; nnnn1 += n1 * (n1 - 1) * n2 * n0;
		rryy2 += r2 * (r2 - 1) * y0 * y1; rrnn2 += r2 * (r2 - 1) * n0 * n1; nnyy2 += n2 * (n2 - 1) * y0 * y1; nnnn2 += n2 * (n2 - 1) * n0 * n1;
	}
	return  (((nnnn0 != 0) ? ((score_t) rrnn0) * ((score_t) nnyy0) / ((score_t) nnnn0) : 0)
		   + ((nnnn1 != 0) ? ((score_t) rrnn1) * ((score_t) nnyy1) / ((score_t) nnnn1) : 0)
		   + ((nnnn2 != 0) ? ((score_t) rrnn2) * ((score_t) nnyy2) / ((score_t) nnnn2) : 0)) - (rryy0 + rryy1 + rryy2) * (1 - 2 * R + 2 * R * R);
}

struct TripartitionInitializer{
	int nThreads = 1;
	vector<score_t> pi;
	vector<vector<vector<int> > > seq;
};

struct Tripartition{
	int nThreads, n = 0;
	vector<vector<array<array<int, 2>, 3> > > cnt; // gene * pos * part * letter
	vector<vector<array<int, 2> > > total; // gene * pos * letter
	const vector<score_t> pi; // gene
	const vector<vector<vector<int> > > seq; // gene * taxon * pos -> letter
	// i, j, k -> taxon; p, q, r -> pos; x, y, z -> part; a, b, c -> letter; g -> gene
	
	Tripartition(const TripartitionInitializer &init): nThreads(init.nThreads), seq(init.seq), pi(init.pi){
		for (const vector<vector<int> > s: seq){
			cnt.emplace_back(s[0].size());
			total.emplace_back(s[0].size());
		}
	}
	
	void resetWork(int start, int end){
		for (int g = start; g < end; g++){
			for (int p = 0; p < total[g].size(); p++){
				for (int a = 0; a < 2; a++){
					cnt[g][p][0][a] = total[g][p][a];
					cnt[g][p][1][a] = 0;
					cnt[g][p][2][a] = 0;
				}
			}
		}
	}
	
	void addTotalWork(int i, int start, int end){
		for (int g = start; g < end; g++){
			for (int p = 0; p < total[g].size(); p++){
				int a = seq[g][i][p];
				if (a != -1) total[g][p][a]++;
			}
		}
	}
	
	void rmvTotalWork(int i, int start, int end){
		for (int g = start; g < end; g++){
			for (int p = 0; p < total[g].size(); p++){
				int a = seq[g][i][p];
				if (a != -1) total[g][p][a]--;
			}
		}
	}
	
	void addWork(int x, int i, int start, int end){
		for (int g = start; g < end; g++){
			for (int p = 0; p < total[g].size(); p++){
				int a = seq[g][i][p];
				if (a != -1) cnt[g][p][x][a]++;
			}
		}
	}
	
	void rmvWork(int x, int i, int start, int end){
		for (int g = start; g < end; g++){
			for (int p = 0; p < total[g].size(); p++){
				int a = seq[g][i][p];
				if (a != -1) cnt[g][p][x][a]--;
			}
		}
	}
	
	void scoreWork(score_t &res, int start, int end){
		score_t temp = 0;
		for (int g = start; g < end; g++){
			temp += scoreGene(cnt[g], pi[g]);
		}
		res = temp;
	}
	
	void reset(){
		vector<thread> thrds;
		for (int t = 1; t < nThreads; t++) thrds.emplace_back(&Tripartition::resetWork, this, pi.size() * t / nThreads, pi.size() * (t + 1) / nThreads);
		resetWork(0, pi.size() / nThreads);
		for (thread &t: thrds) t.join();
	}
	
	void addTotal(int i){
		n++;
		vector<thread> thrds;
		for (int t = 1; t < nThreads; t++) thrds.emplace_back(&Tripartition::addTotalWork, this, i, pi.size() * t / nThreads, pi.size() * (t + 1) / nThreads);
		addTotalWork(i, 0, pi.size() / nThreads);
		for (thread &t: thrds) t.join();
	}
	
	void rmvTotal(int i){
		n--;
		vector<thread> thrds;
		for (int t = 1; t < nThreads; t++) thrds.emplace_back(&Tripartition::rmvTotalWork, this, i, pi.size() * t / nThreads, pi.size() * (t + 1) / nThreads);
		rmvTotalWork(i, 0, pi.size() / nThreads);
		for (thread &t: thrds) t.join();
	}
	
	void add(int x, int i){
		vector<thread> thrds;
		for (int t = 1; t < nThreads; t++) thrds.emplace_back(&Tripartition::addWork, this, x, i, pi.size() * t / nThreads, pi.size() * (t + 1) / nThreads);
		addWork(x, i, 0, pi.size() / nThreads);
		for (thread &t: thrds) t.join();
	}
	
	void rmv(int x, int i){
		vector<thread> thrds;
		for (int t = 1; t < nThreads; t++) thrds.emplace_back(&Tripartition::rmvWork, this, x, i, pi.size() * t / nThreads, pi.size() * (t + 1) / nThreads);
		rmvWork(x, i, 0, pi.size() / nThreads);
		for (thread &t: thrds) t.join();
	}
	
	score_t score(){
		vector<score_t> res(nThreads);
		vector<thread> thrds;
		for (int t = 1; t < nThreads; t++) thrds.emplace_back(&Tripartition::scoreWork, this, ref(res[t]), pi.size() * t / nThreads, pi.size() * (t + 1) / nThreads);
		scoreWork(res[0], 0, pi.size() / nThreads);
		for (thread &t: thrds) t.join();
		for (int t = 1; t < nThreads; t++) res[0] += res[t];
		return res[0];
	}
};

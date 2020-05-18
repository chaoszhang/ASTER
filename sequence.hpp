typedef __float128 score_t;

#include<vector>
#include<array>

using namespace std;

inline score_t scorePos(const array<array<int, 4>, 3> cnt, const array<score_t, 4> pi){
	int sum[3] = {0, 0, 0};
	score_t res = 0;
	for (int x = 0; x < 3; x++){
		for (int a = 0; a < 4; a++){
			sum[x] += cnt[x][a];
		}
	}
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
}

struct Tripartition{
	int npos, ntaxa;
	vector<array<array<int, 4>, 3> > cnt; // pos * part * letter
	vector<array<int, 4> > total; // pos * letter
	vector<array<score_t, 4> > pi; // pos * letter
	vector<vector<int> > seq; // taxon * pos -> letter
	// i, j, k -> taxon; p, q, r -> pos; x, y, z -> part; a, b, c -> letter
	
	void reset(){
		for (int p = 0; p < npos; p++){
			for (int a = 0; a < 4; a++){
				cnt[p][0][a] = total[p][a];
				cnt[p][1][a] = 0;
				cnt[p][2][a] = 0;
			}
		}
	}
	
	void addTotal(int i){
		for (int p = 0; p < npos; p++){
			int a = seq[i][p];
			if (a != -1) total[p][a]++;
		}
	}
	
	void add(int x, int i){
		for (int p = 0; p < npos; p++){
			int a = seq[i][p];
			if (a != -1) cnt[p][x][a]++;
		}
	}
	
	void rmv(int x, int i){
		for (int p = 0; p < npos; p++){
			int a = seq[i][p];
			if (a != -1) cnt[p][x][a]--;
		}
	}
	
	score_t score(){
		score_t res = 0;
		for (int p = 0; p < npos; p++){
			res += scorePos(cnt[p], pi[p]);
		}
		return res;
	}
};

#define OBJECTIVE_VERSION "0"

#include<vector>
#include<cmath>
#include "threadpool.hpp"

using namespace std;

struct TripartitionInitializer{
    vector<vector<double> > M, cnt;
};

struct Tripartition{
	vector<char> color;
	vector<vector<double> > &M, &cnt;
    vector<double> S[3], C[3], Q[3];
    int n, m;


	Tripartition(TripartitionInitializer &init): M(init.M), cnt(init.cnt), n(init.M.size()), m(init.M[0].size()), color(init.M.size(), -1){
        for (int x = 0; x < 3; x++){
            S[x].resize(m);
            C[x].resize(m);
            Q[x].resize(m);
        }
	}

	void updatePart(int part, int x, int i){
        if (part) return;
		int y = color[i];
		color[i] = x;
		for (int j = 0; j < m; j++){
            if (y != -1){
                C[y][j] -= cnt[i][j];
                Q[y][j] -= cnt[i][j] * cnt[i][j];
                S[y][j] -= M[i][j];
            }
            if (x != -1){
                C[x][j] += cnt[i][j];
                Q[x][j] += cnt[i][j] * cnt[i][j];
                S[x][j] += M[i][j];
            }
        }
	}
	
	score_t scorePart(int part){
		if (part) return 0;
        double res = 0;
        for (int j = 0; j < m; j++){
            res += (C[0][j] * C[0][j] - Q[0][j]) * S[1][j] * S[2][j]
                 + (C[1][j] * C[1][j] - Q[1][j]) * S[2][j] * S[0][j]
                 + (C[2][j] * C[2][j] - Q[2][j]) * S[0][j] * S[1][j];
        }
        return res;
	}
};
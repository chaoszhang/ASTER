#define OBJECTIVE_VERSION "0"

#include<vector>
#include<array>
#include<cmath>
#include "threadpool.hpp"

#define CUSTOMIZED_ANNOTATION
#define LOCAL_BOOTSTRAP

using namespace std;

struct CustomizedAnnotation{
	vector<array<score_t,3> > bs;

	CustomizedAnnotation(){}
	CustomizedAnnotation(int len): bs(len){}

	array<int, 3> bootstrap(){
		array<int, 3> result = {};
		for (const array<score_t,3> &e: bs){
			if (e[0] > e[1] && e[0] > e[2]) result[0]++;
			if (e[1] > e[0] && e[1] > e[2]) result[1]++;
			if (e[2] > e[0] && e[2] > e[1]) result[2]++;
		}
		return result;
	}

	CustomizedAnnotation operator+ (const CustomizedAnnotation& o) const{
		CustomizedAnnotation r(bs.size());
		for (int i = 0; i < bs.size(); i++){
			for (int j = 0; j < 3; j++){
				r.bs[i][j] = bs[i][j] + o.bs[i][j];
			}
		}
		return r;
	}

	CustomizedAnnotation& operator+= (const CustomizedAnnotation& o){
		for (int i = 0; i < bs.size(); i++){
			for (int j = 0; j < 3; j++){
				bs[i][j] += o.bs[i][j];
			}
		}
		return *this;
	}
};

struct TripartitionInitializer{
    vector<vector<double> > M, cnt;
    int nThreads = 1;
};

struct Tripartition{
	vector<vector<char> > color;
	vector<vector<double> > &M, &cnt;
    vector<double> S[3], C[3], Q[3];
    int n, m, nThreads;

	Tripartition(TripartitionInitializer &init): M(init.M), cnt(init.cnt), n(init.M.size()), m(init.M[0].size()), nThreads(init.nThreads), color(init.nThreads, vector<char>(init.M.size(), -1)){
        for (int x = 0; x < 3; x++){
            S[x].resize(m);
            C[x].resize(m);
            Q[x].resize(m);
        }
	}

	void updatePart(int part, int x, int i){
        int y = color[part][i];
		color[part][i] = x;
		for (int j = m * part / nThreads; j < m * (part + 1) / nThreads; j++){
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
		double res = 0;
        for (int j = m * part / nThreads; j < m * (part + 1) / nThreads; j++){
            res += (C[0][j] * C[0][j] - Q[0][j]) * S[1][j] * S[2][j]
                 + (C[1][j] * C[1][j] - Q[1][j]) * S[2][j] * S[0][j]
                 + (C[2][j] * C[2][j] - Q[2][j]) * S[0][j] * S[1][j];
        }
        return res;
	}
};

struct Quadrupartition{
    const int ufb_size = 1000;
	vector<char> color;
	vector<vector<double> > &M, &cnt;
    vector<double> S[4], C[4];
    vector<vector<int> > bsCol;
    int n, m, nThreads;

	Quadrupartition(TripartitionInitializer &init): M(init.M), cnt(init.cnt), n(init.M.size()), m(init.M[0].size()), nThreads(init.nThreads), color(init.M.size(), -1), bsCol(ufb_size){
        for (int x = 0; x < 4; x++){
            S[x].resize(m);
            C[x].resize(m);
        }
        for (int i = 0; i < ufb_size; i++){
            for (int j = 0; j < m; j++){
                bsCol[i].push_back((rand() * 32768LL + rand()) % m);
            }
        }
	}

	void update(int x, int i){
        int y = color[i];
		color[i] = x;
		for (int j = 0; j < m; j++){
            if (y != -1){
                C[y][j] -= cnt[i][j];
                S[y][j] -= M[i][j];
            }
            if (x != -1){
                C[x][j] += cnt[i][j];
                S[x][j] += M[i][j];
            }
        }
	}

    array<double, 3> score(){
        array<double, 3> res;
		vector<double> t1(m), t2(m), t3(m);
        for (int j = 0; j < m; j++){
            res[0] += C[0][j] * C[1][j] * S[2][j] * S[3][j] + S[0][j] * S[1][j] * C[2][j] * C[3][j];
            res[1] += C[0][j] * C[2][j] * S[1][j] * S[3][j] + S[0][j] * S[2][j] * C[1][j] * C[3][j];
            res[2] += C[0][j] * C[3][j] * S[1][j] * S[2][j] + S[0][j] * S[3][j] * C[1][j] * C[2][j];
        }
        return res;
	}
	
	CustomizedAnnotation annotate(){
        CustomizedAnnotation res(ufb_size);
		vector<double> t1(m), t2(m), t3(m);
        for (int j = 0; j < m; j++){
            t1[j] = C[0][j] * C[1][j] * S[2][j] * S[3][j] + S[0][j] * S[1][j] * C[2][j] * C[3][j];
            t2[j] = C[0][j] * C[2][j] * S[1][j] * S[3][j] + S[0][j] * S[2][j] * C[1][j] * C[3][j];
            t3[j] = C[0][j] * C[3][j] * S[1][j] * S[2][j] + S[0][j] * S[3][j] * C[1][j] * C[2][j];
        }
        for (int i = 0; i < ufb_size; i++){
            for (int j: bsCol[i]){
                res.bs[i][0] += t1[j];
                res.bs[i][1] += t2[j];
                res.bs[i][2] += t3[j];
            }
        }
        return res;
	}
};
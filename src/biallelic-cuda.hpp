#define OBJECTIVE_VERSION "0"

#include<iostream>
#include<vector>
#include<array>
#include<algorithm>
#include<bitset>
//#include "threadpool.hpp"
#include "biallelic-cuda.cuh"

using namespace std;

#define CACHE_ALIGNMENT 64

struct TripartitionInitializer{
	vector<int> breakPoints;
	vector<array<vector<bool>, 4> > h;
	vector<array<score_t, 4> > p;
	vector<score_t> w;
	int numThreads = 1;
};

struct Tripartition{
	GroupIBatch *GI;
	GroupFBatch *GF;
	SiteBatch *S;
	HotBatch *H;
	double *BS;
	int n, kc;
	vector<char> color;
	vector<double> batchScore;

	Tripartition(const TripartitionInitializer &TI): n(TI.h.size()), color(TI.h.size(), -1){
		vector<pair<int, int> > lengthSorted;
		int k = TI.w.size();
		kc = (k+BATCH-1)/BATCH;
		for (int i = 0; i < k; i++){
			lengthSorted.push_back({i, TI.breakPoints[i] - (int)((i == 0) ? 0 : TI.breakPoints[i - 1])});
		}
		sort(lengthSorted.begin(), lengthSorted.end(), [](const pair<int, int> &a, const pair<int, int> &b)->bool{return a.second < b.second;});
		vector<array<int,2> > VGI;
		vector<array<array<double,BATCH>,5> > VGF;
		vector<BATCH_TYPE> VH;
		int vsSize = 0;
		for (int t = 0; t < kc; t++){
			int maxNSites = 0, id[BATCH] = {};
			VGI.emplace_back();
			VGF.emplace_back();
			VGI.back()[0] = vsSize;
			for (int b = 0; b < BATCH; b++){
				if (t * BATCH + b >= k) break;
				int i = lengthSorted[t * BATCH + b].first;
				id[b] = i;
				int start = ((i == 0) ? 0 : TI.breakPoints[i - 1]);
				int end = TI.breakPoints[i];
				for (int p = 0; p < 4; p++) VGF.back()[p][b] = TI.p[i][p];
				VGF.back()[4][b] = TI.w[i];
				maxNSites = max(maxNSites, end - start);
			}
			vsSize += maxNSites;
			VGI.back()[1] = vsSize;
			for (int r = 0; r < n; r++){
				for (int j = 0; j < maxNSites; j++){
					for (int p = 0; p < 4; p++){
						VH.push_back(0);
						for (int b = 0; b < BATCH; b++){
							int i = id[b];
							int start = ((i == 0) ? 0 : TI.breakPoints[i - 1]);
							int end = TI.breakPoints[i];
							if (start + j >= end) continue;
							VH.back() |= (TI.h[r][p][start + j] << b);
						}
					}
				}
			}
		}
		vector<SiteBatch> VS(vsSize);
		batchScore.resize(kc);
		myCudaMalloc((void**) &GI, VGI.size() * sizeof(GroupIBatch));
		myCudaMalloc((void**) &GF, VGF.size() * sizeof(GroupFBatch));
		myCudaMalloc((void**) &S, VS.size() * sizeof(SiteBatch));
		myCudaMalloc((void**) &H, VH.size() * sizeof(BATCH_TYPE));
		myCudaMalloc((void**) &BS, batchScore.size() * sizeof(double));
		myCudaMemcpyH2D(GI, &VGI[0], VGI.size() * sizeof(GroupIBatch));
		myCudaMemcpyH2D(GF, &VGF[0], VGF.size() * sizeof(GroupFBatch));
		myCudaMemcpyH2D(S, &VS[0], VS.size() * sizeof(SiteBatch));
		myCudaMemcpyH2D(H, &VH[0], VH.size() * sizeof(BATCH_TYPE));
		myCudaMemcpyH2D(BS, &batchScore[0], batchScore.size() * sizeof(double));
		//for (auto e: VH) printf("%x ",e);
	}

	~Tripartition(){
		myCudaFree(GI);
		myCudaFree(GF);
		myCudaFree(S);
		myCudaFree(H);
		myCudaFree(BS);
	}

	void update(int x, int i){
		int y = color[i];
		cudaUpdate(GI, S, H, n, x, y, i, kc, BATCH);
		SiteBatch s;
		myCudaMemcpyD2H(&s, S, sizeof(SiteBatch));
		//cerr << x << " " << y << endl;
		//cerr << s[0][0][0] << " " << s[1][0][0] << " " << s[2][0][0] << " " << s[3][0][0] << endl;
		//cerr << s[0][1][0] << " " << s[1][1][0] << " " << s[2][1][0] << " " << s[3][1][0] << endl;
		//cerr << s[0][2][0] << " " << s[1][2][0] << " " << s[2][2][0] << " " << s[3][2][0] << endl << endl;
		color[i] = x;
	}

	score_t score(){
		cudaScore(GI, GF, S, BS, kc, BATCH);
		myCudaMemcpyD2H(&batchScore[0], BS, kc * sizeof(double));
		score_t result = 0.0;
		for (double s: batchScore) result += s;
		//cerr << result << endl;
		return result;
	}
};
#define OBJECTIVE_VERSION "0"

#include<vector>
#include<array>
#include<thread>
#include "threadpool.hpp"

using namespace std;

inline long long scorePair(long long oddX0, long long oddX1, long long oddX2, long long oddY0, long long oddY1, long long oddY2, score_t piXpiYodd,
		long long evenX0, long long evenX1, long long evenX2, long long evenY0, long long evenY1, long long evenY2, score_t piXpiYeven){
	long long odd0 = oddX0 + oddY0, odd1 = oddX1 + oddY1, odd2 = oddX2 + oddY2;
	long long even0 = evenX0 + evenY0, even1 = evenX1 + evenY1, even2 = evenX2 + evenY2;
	return (piXpiYodd * odd0 * (odd0 - 1) - oddX0 * oddY0) * (2 * piXpiYeven * even1 * even2 - evenX1 * evenY2 - evenY1 * evenX2)
		 + (piXpiYeven * even0 * (even0 - 1) - evenX0 * evenY0) * (2 * piXpiYodd * odd1 * odd2 - oddX1 * oddY2 - oddY1 * oddX2);
}

inline score_t scorePos(const array<int, 3> oddX, const array<int, 3> oddY, score_t piXpiYodd,
		const array<int, 3> evenX, const array<int, 3> evenY, score_t piXpiYeven){
	return scorePair(oddX[0], oddX[1], oddX[2], oddY[0], oddY[1], oddY[2], piXpiYodd,
					 evenX[0], evenX[1], evenX[2], evenY[0], evenY[1], evenY[2], piXpiYeven)
		 + scorePair(oddX[1], oddX[2], oddX[0], oddY[1], oddY[2], oddY[0], piXpiYodd,
					 evenX[1], evenX[2], evenX[0], evenY[1], evenY[2], evenY[0], piXpiYeven)
		 + scorePair(oddX[2], oddX[0], oddX[1], oddY[2], oddY[0], oddY[1], piXpiYodd,
					 evenX[2], evenX[0], evenX[1], evenY[2], evenY[0], evenY[1], piXpiYeven);
}

struct TripartitionInitializer{
	struct Tree{
		int npos = 0, ntree = 0;
		vector<score_t> piXpiY;
		vector<vector<bool> > seqX, seqY;
		vector<vector<int> > speciesMapper;

		Tree(){}
		Tree(int n): seqX(n), seqY(n), speciesMapper(n){}
	};
	int nThreads = 1;
	vector<vector<Tree> > trees;
};

struct Tripartition{
	struct Tree{
		vector<vector<array<int, 3> > > cntX, cntY; // tree * pos * part
		vector<score_t> scores; // tree;
		vector<bool> valid; // tree;

		Tree(){}
		Tree(int ntree, int npos): scores(ntree), valid(ntree),
			cntX(ntree, vector<array<int, 3> >(npos)), cntY(ntree, vector<array<int, 3> >(npos)){}
	};

	const TripartitionInitializer& TI;
	vector<vector<int> > color;
	vector<vector<Tree> > trees; // gene * bin

	Tripartition(const TripartitionInitializer &init): TI(init), trees(init.trees.size()),
			color(init.nThreads, vector<int>(init.trees[0][0].seqX.size(), -1)){
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
					for (int p = 0; p < trees[a][b].cntX[t].size(); p++){
						if (y != -1) {
							trees[a][b].cntX[t][p][y] -= TI.trees[a][b].seqX[i][p];
							trees[a][b].cntY[t][p][y] -= TI.trees[a][b].seqY[i][p];
						}
						if (x != -1) {
							trees[a][b].cntX[t][p][x] += TI.trees[a][b].seqX[i][p];
							trees[a][b].cntY[t][p][x] += TI.trees[a][b].seqY[i][p];
						}
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
					for (int p = 0; p < trees[a][b].cntX[t].size(); p += 2){
						temp += scorePos(trees[a][b].cntX[t][p], trees[a][b].cntY[t][p], TI.trees[a][b].piXpiY[p],
										 trees[a][b].cntX[t][p+1], trees[a][b].cntY[t][p+1], TI.trees[a][b].piXpiY[p+1]);
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

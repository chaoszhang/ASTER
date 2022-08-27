#define OBJECTIVE_VERSION "0"

#include<vector>
#include<array>
#include<thread>
#include "threadpool.hpp"

#define G_SUPPORT
#define CUSTOMIZED_LENGTH
using namespace std;

inline score_t scorePair(long long oddX0, long long oddX1, long long oddX2, long long oddY0, long long oddY1, long long oddY2, score_t piXpiYodd,
		long long evenX0, long long evenX1, long long evenX2, long long evenY0, long long evenY1, long long evenY2, score_t piXpiYeven){
	long long odd0 = oddX0 + oddY0, odd1 = oddX1 + oddY1, odd2 = oddX2 + oddY2;
	long long even0 = evenX0 + evenY0, even1 = evenX1 + evenY1, even2 = evenX2 + evenY2;
	return (piXpiYodd * odd0 * (odd0 - 1) - oddX0 * oddY0) * (2 * piXpiYeven * even1 * even2 - evenX1 * evenY2 - evenY1 * evenX2)
		 + (piXpiYeven * even0 * (even0 - 1) - evenX0 * evenY0) * (2 * piXpiYodd * odd1 * odd2 - oddX1 * oddY2 - oddY1 * oddX2);
}

inline score_t scorePos(const array<unsigned short, 3> oddX, const array<unsigned short, 3> oddY, score_t piXpiYodd,
		const array<unsigned short, 3> evenX, const array<unsigned short, 3> evenY, score_t piXpiYeven){
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
		vector<float> piXpiY;
		vector<vector<bool> > seqX, seqY;
		vector<vector<unsigned short> > speciesMapper;

		Tree(){}
		Tree(int n): seqX(n), seqY(n), speciesMapper(n){}
	};
	int nThreads = 1;
	vector<vector<Tree> > trees;
};

struct Tripartition{
	struct Tree{
		vector<vector<array<unsigned short, 3> > > cntX, cntY; // tree * pos * part
		vector<score_t> scores; // tree;
		vector<bool> valid; // tree;

		Tree(){}
		Tree(int ntree, int npos): scores(ntree), valid(ntree),
			cntX(ntree, vector<array<unsigned short, 3> >(npos)), cntY(ntree, vector<array<unsigned short, 3> >(npos)){}
	};

	const TripartitionInitializer& TI;
	vector<vector<char> > color;
	vector<vector<Tree> > trees; // gene * bin

	Tripartition(const TripartitionInitializer &init): TI(init), trees(init.trees.size()),
			color(init.nThreads, vector<char>(init.trees[0][0].seqX.size(), -1)){
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

inline score_t supportPair(long long oddX0, long long oddX1, long long oddX2, long long oddX3, 
		long long oddY0, long long oddY1, long long oddY2, long long oddY3, score_t piXpiYodd,
		long long evenX0, long long evenX1, long long evenX2, long long evenX3,
		long long evenY0, long long evenY1, long long evenY2, long long evenY3, score_t piXpiYeven){
	long long odd0 = oddX0 + oddY0, odd1 = oddX1 + oddY1, odd2 = oddX2 + oddY2, odd3 = oddX3 + oddY3;
	long long even0 = evenX0 + evenY0, even1 = evenX1 + evenY1, even2 = evenX2 + evenY2, even3 = evenX3 + evenY3;
	return (2 * piXpiYodd * odd0 * odd1 - oddX0 * oddY1 - oddY0 * oddX1) * (2 * piXpiYeven * even2 * even3 - evenX2 * evenY3 - evenY2 * evenX3)
		 + (2 * piXpiYeven * even0 * even1 - evenX0 * evenY1 - evenY0 * evenX1) * (2 * piXpiYodd * odd2 * odd3 - oddX2 * oddY3 - oddY2 * oddX3);
}

array<double, 3> supportPos(const array<unsigned short, 4> oddX, const array<unsigned short, 4> oddY, score_t piXpiYodd,
		const array<unsigned short, 4> evenX, const array<unsigned short, 4> evenY, score_t piXpiYeven){
	return {(double) supportPair(oddX[0], oddX[1], oddX[2], oddX[3], oddY[0], oddY[1], oddY[2], oddY[3], piXpiYodd,
					    		 evenX[0], evenX[1], evenX[2], evenX[3], evenY[0], evenY[1], evenY[2], evenY[3], piXpiYeven),
		    (double) supportPair(oddX[0], oddX[2], oddX[1], oddX[3], oddY[0], oddY[2], oddY[1], oddY[3], piXpiYodd,
					             evenX[0], evenX[2], evenX[1], evenX[3], evenY[0], evenY[2], evenY[1], evenY[3], piXpiYeven),
		    (double) supportPair(oddX[0], oddX[3], oddX[1], oddX[2], oddY[0], oddY[3], oddY[1], oddY[2], piXpiYodd,
					             evenX[0], evenX[3], evenX[1], evenX[2], evenY[0], evenY[3], evenY[1], evenY[2], piXpiYeven)};
}

struct Quadrupartition{
	struct Tree{
		vector<vector<array<unsigned short, 4> > > cntX, cntY; // tree * pos * part
		vector<array<score_t, 3> > scores; // tree;
		vector<bool> valid; // tree;

		Tree(){}
		Tree(int ntree, int npos): scores(ntree), valid(ntree),
			cntX(ntree, vector<array<unsigned short, 4> >(npos)), cntY(ntree, vector<array<unsigned short, 4> >(npos)){}
	};

	const TripartitionInitializer& TI;
	vector<vector<char> > color;
	vector<vector<Tree> > trees; // gene * bin

	Quadrupartition(const TripartitionInitializer &init): TI(init), trees(init.trees.size()),
			color(init.nThreads, vector<char>(init.trees[0][0].seqX.size(), -1)){
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
	
	array<double, 3> scorePart(int part){
		int start = trees.size() * part / TI.nThreads, end = trees.size() * (1 + part) / TI.nThreads;
		array<double, 3> result;
		result[0] = 0;
		result[1] = 0;
		result[2] = 0;
		for (int a = start; a < end; a++){
			for (int b = 0; b < trees[a].size(); b++){
				for (int t = 0; t < trees[a][b].scores.size(); t++){
					if (trees[a][b].valid[t]){
						result[0] += trees[a][b].scores[t][0];
						result[1] += trees[a][b].scores[t][1];
						result[2] += trees[a][b].scores[t][2];
						continue;
					}
					score_t temp0 = 0, temp1 = 0, temp2 = 0;
					for (int p = 0; p < trees[a][b].cntX[t].size(); p += 2){
						array<double, 3> cnt = supportPos(trees[a][b].cntX[t][p], trees[a][b].cntY[t][p], TI.trees[a][b].piXpiY[p],
										                  trees[a][b].cntX[t][p+1], trees[a][b].cntY[t][p+1], TI.trees[a][b].piXpiY[p+1]);
						temp0 += cnt[0];
						temp1 += cnt[1];
						temp2 += cnt[2];
					}
					trees[a][b].scores[t][0] = temp0;
					trees[a][b].scores[t][1] = temp1;
					trees[a][b].scores[t][2] = temp2;
					trees[a][b].valid[t] = true;
					result[0] += trees[a][b].scores[t][0];
					result[1] += trees[a][b].scores[t][1];
					result[2] += trees[a][b].scores[t][2];
				}
			}
		}
		return result;
	}

	void update(int x, int i){
		vector<thread> thrds;
		for (int p = 1; p < TI.nThreads; p++) thrds.emplace_back(&Quadrupartition::updatePart, this, p, x, i);
		updatePart(0, x, i);
		for (thread &t: thrds) t.join();
	}
	
	array<double, 3> score(){
		vector<thread> thrds;
		vector<future<array<double, 3> > > temp(TI.nThreads);
		for (int p = 1; p < TI.nThreads; p++){
			temp[p] = async(launch::async, &Quadrupartition::scorePart, this, p);
		}
		array<double, 3> res = scorePart(0);
		for (int p = 1; p < TI.nThreads; p++){
			array<double, 3> t = temp[p].get();
			res[0] += t[0];
			res[1] += t[1];
			res[2] += t[2];
		}
		return res;
	}

	string length(const array<double, 3> &scores){
		if (scores[0] <= 0) return to_string(-0.0);
		if (scores[0] < max(scores[1], scores[2])) return to_string(-0.0);
		if (max(scores[1], scores[2]) <= 0) return "Inf";
		return to_string((log(scores[0]) - log(max(scores[1], scores[2]))) / 2);
	}
};

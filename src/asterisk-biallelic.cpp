#define DRIVER_VERSION "1"

#include<iostream>
#include<fstream>
#include<unordered_map>
#include<unordered_set>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<algorithm>
#include<random>
#include<thread>
#include<mutex>

#define LARGE_DATA
#ifdef LARGE_DATA
typedef long double score_t;
typedef long long count_t;
#else
typedef double score_t;
typedef int count_t;
#endif

#define ERROR_TOLERANCE 0.001

#include "argparser.hpp"
//#define USE_CUDA
#ifdef USE_CUDA
#include "biallelic-cuda.hpp"
#else
#include "biallelic.hpp"
#endif

#include "algorithms.hpp"

using namespace std;

const int NUM_BINS = 20, MIN_JURY_SIZE = 50, MAX_CANDIDATE_SIZE = 200, VIP_CANDIDATE_SIZE = 100, MIN_BIN_SIZE = 10;
const double KEEP_RATE = 1;

MetaAlgorithm meta;
TripartitionInitializer &tripInit = meta.tripInit;

vector<string> &names = meta.names;
unordered_map<string, int> &name2id = meta.name2id;

mutex mtx;

string formatName(const string name){
	string res;
	for (char c: name){
		if (c != '>' && c != ' ' && c != '\t') res += c;
	}
	return res;
}

vector<pair<int, double> > infoSort(vector<array<vector<bool>, 4> > &sites, vector<int> &vips){
	mt19937 generator;
	vector<int> candidate, ordering, jury;
	unordered_set<int> vipset;
	for (int i: vips){
		ordering.push_back(i);
	}
	std::shuffle(ordering.begin(), ordering.end(), generator);
	for (int i = 0, imax = min(VIP_CANDIDATE_SIZE, (int) vips.size()); i < imax; i++){
		candidate.push_back(ordering[i]);
		vipset.insert(i);
	}
	ordering.clear();
	for (int i = 0; i < sites.size(); i++){
		if (vipset.count(i) == 0) ordering.push_back(i);
	}
	shuffle(ordering.begin(), ordering.end(), generator);
	for (int i = 0, imax = min((int) ordering.size(), MAX_CANDIDATE_SIZE - (int) candidate.size()); i < imax; i++){
		candidate.push_back(ordering[i]);
	}
	vector<pair<int, double> > candidatePair, allPair;
	for (int i: candidate){
		double MI = 0;
		for (int j: candidate){
			if (i == j) continue;
			double freq[4][4] = {}, ifreq[4] = {}, jfreq[4] = {}, tfreq = 0;
			for (int x = 0; x < 4; x++){
				for (int y = 0; y < 4; y++){
					int cnt = 0;
					for (int p = 0; p < sites[i][x].size(); p++){
						cnt += (sites[i][x][p] & sites[j][y][p]);
					}
					freq[x][y] = cnt;
					ifreq[x] += freq[x][y];
					jfreq[y] += freq[x][y]; 
					tfreq += freq[x][y];
				}
			}
			for (int x = 0; x < 4; x++){
				for (int y = 0; y < 4; y++){
					if (freq[x][y] == 0) continue;
					MI += freq[x][y] * (log(freq[x][y] * tfreq / (ifreq[x] * jfreq[y])));
				}
			}
		}
		candidatePair.push_back({i, MI});
	}
	sort(candidatePair.begin(), candidatePair.end(), [](pair<int, double> a, pair<int, double> b){return a.second > b.second;});
	int jurySize = max(MIN_JURY_SIZE, MAX_CANDIDATE_SIZE * MAX_CANDIDATE_SIZE / (int) sites.size());
	for (int i = 0; i < min((int) candidate.size(), jurySize); i++){
		jury.push_back(candidatePair[i].first);
	}
	for (int i = 0; i < sites.size(); i++){
		double MI = 0;
		int juryCnt = 0;
		for (int j: jury){
			if (i == j) continue;
			juryCnt++;
			double freq[4][4] = {}, ifreq[4] = {}, jfreq[4] = {}, tfreq = 0;
			for (int x = 0; x < 4; x++){
				for (int y = 0; y < 4; y++){
					int cnt = 0;
					for (int p = 0; p < sites[i][x].size(); p++){
						cnt += (sites[i][x][p] & sites[j][y][p]);
					}
					freq[x][y] = cnt;
					ifreq[x] += freq[x][y];
					jfreq[y] += freq[x][y];
					tfreq += freq[x][y];
				}
			}
			for (int x = 0; x < 4; x++){
				for (int y = 0; y < 4; y++){
					if (freq[x][y] == 0) continue;
					MI += freq[x][y] / tfreq * log(freq[x][y] * tfreq / (ifreq[x] * jfreq[y]));
				}
			}
		}
		allPair.push_back({i, MI / juryCnt});
	}
	sort(allPair.begin(), allPair.end(), [](pair<int, double> a, pair<int, double> b){return a.second > b.second;});
	//for (auto p:allPair) cerr << p.first << ":" << p.second << "\t";
	return allPair;
}

void addBin(const vector<int> &bin, const array<score_t, 4> &PI, const vector<string> &geneNames, const vector<string> &geneSeqs){
	if (bin.size() < 2) return;
	const lock_guard<mutex> lock(mtx);
	int oldL = (tripInit.breakPoints.size()) ? tripInit.breakPoints.back() : 0;
	int L = oldL + bin.size();
	tripInit.breakPoints.push_back(L);
	tripInit.p.push_back(PI);
	for (int I = 0; I < tripInit.h.size(); I++){
		for (int x = 0; x < 4; x++){
			tripInit.h[I][x].resize(L, false);
		}
	}
	for (int i = 0; i < geneSeqs.size(); i++){
		int I = name2id[geneNames[i]];
		for (int j = 0; j < bin.size(); j++){
			switch (geneSeqs[i][bin[j]]){
				case 'A': case 'a': tripInit.h[I][0][oldL + j] = true; break; 
				case 'C': case 'c': tripInit.h[I][1][oldL + j] = true; break;
				case 'G': case 'g': tripInit.h[I][2][oldL + j] = true; break;
				case 'T': case 't': tripInit.h[I][3][oldL + j] = true; break;
			}
		}
	}
	tripInit.w.push_back(1.0 / (bin.size() - 1));
}

void formatGene(const vector<string> geneNames, const vector<string> geneSeqs){
	{
		const lock_guard<mutex> lock(mtx);
		int L = (tripInit.breakPoints.size()) ? tripInit.breakPoints.back() : 0;
		for (string n: geneNames){
			if (name2id.count(n) == 0){
				int I = names.size();
				name2id[n] = I;
				names.push_back(n);
				tripInit.h.emplace_back();
				for (int x = 0; x < 4; x++){
					tripInit.h[I][x].resize(L, false);
				}
			}
		}
	}
	vector<int> triallelicSingle, triallelicNoSingle;
	/*
	vector<int> biallelic, triallelic;
	vector<int> bitriallelic, quadriallelic;
	vector<int> quadriallelicSingle, quadriallelicNoSingle;
	vector<int> quadriallelicFew, quadriallelicMany;
	vector<int> quadriallelicNotRapid, quadriallelicRapid;
	*/
	unordered_set<int> biallelicSet, triallelicSet, quadriallelicSet;
	unordered_set<int> bitriallelicSet;
	int totalA = 0, totalC = 0, totalG = 0, totalT = 0;
	vector<array<vector<bool>, 4> > sites(geneSeqs[0].size());
	for (int j = 0; j < geneSeqs[0].size(); j++){
		int cntA = 0, cntC = 0, cntG = 0, cntT = 0;
		for (int x = 0; x < 4; x++){
			sites[j][x].resize(geneSeqs.size(), false);
		}
		for (int i = 0; i < geneSeqs.size(); i++){
			switch (geneSeqs[i][j]){
				case 'A': case 'a': cntA++; sites[j][0][i] = true; break;
				case 'C': case 'c': cntC++; sites[j][1][i] = true; break;
				case 'G': case 'g': cntG++; sites[j][2][i] = true; break;
				case 'T': case 't': cntT++; sites[j][3][i] = true; break;
			}
		}
		int allelicCnt = (cntA > 0) + (cntC > 0) + (cntG > 0) + (cntT > 0);
		int singleCnt = (cntA == 1) + (cntC == 1) + (cntG == 1) + (cntT == 1);
		int fewCnt = (cntA <= 2) + (cntC <= 2) + (cntG <= 2) + (cntT <= 2);
		int rTh = max(4, (int) geneSeqs.size() / 20);
		bool rapid = (cntA >= rTh && cntC >= rTh && cntG >= rTh && cntT >= rTh);
		if (allelicCnt == 3 && singleCnt > 0) triallelicNoSingle.push_back(j);
		/*
		if (allelicCnt == 2 && singleCnt == 0) biallelic.push_back(j);
		if (allelicCnt == 3) triallelic.push_back(j);
		if (allelicCnt == 3 && singleCnt == 0) triallelicSingle.push_back(j);
		if (allelicCnt == 2 || allelicCnt == 3) bitriallelic.push_back(j);
		if (allelicCnt == 4) quadriallelic.push_back(j);
		if (allelicCnt == 4 && singleCnt > 0) quadriallelicSingle.push_back(j);
		if (allelicCnt == 4 && singleCnt == 0) quadriallelicNoSingle.push_back(j);
		if (allelicCnt == 4 && singleCnt == 0 && fewCnt > 0) quadriallelicFew.push_back(j);
		if (allelicCnt == 4 && fewCnt == 0) quadriallelicMany.push_back(j);
		if (allelicCnt == 4 && fewCnt == 0 && !rapid) quadriallelicNotRapid.push_back(j);
		if (allelicCnt == 4 && rapid) quadriallelicRapid.push_back(j);
		*/
		if (allelicCnt == 2 && singleCnt == 0) biallelicSet.insert(j);
		//if ((allelicCnt == 2 && singleCnt == 0) || allelicCnt == 3) bitriallelicSet.insert(j);
		if (allelicCnt == 3) triallelicSet.insert(j);
		if (allelicCnt == 4) quadriallelicSet.insert(j);
		totalA += cntA; totalC += cntC; totalG += cntG; totalT += cntT; 
	}
	score_t total = totalA + totalC + totalG + totalT;
	array<score_t, 4> PI = {totalA / total, totalC / total, totalG / total, totalT / total};
	vector<pair<int, double> > ordered = infoSort(sites, triallelicNoSingle);
	int numKeep = KEEP_RATE * ordered.size();
	vector<vector<int> > orderedDivision(3);
	for (int i = 0; i < numKeep; i++){
		int p = ordered[i].first;
		if (biallelicSet.count(p)) orderedDivision[0].push_back(p);
		if (triallelicSet.count(p)) orderedDivision[1].push_back(p);
		if (quadriallelicSet.count(p)) orderedDivision[2].push_back(p);
	}
	for (vector<int> &division: orderedDivision){
		int nBin = max(1, min(NUM_BINS, (int) division.size() / MIN_BIN_SIZE));
		for (int i = 0; i < nBin; i++){
			vector<int> bin;
			for (int j = i * division.size() / nBin, jmax = (i + 1) * division.size() / nBin; j < jmax; j++){
				bin.push_back(division[j]);
			}
			addBin(bin, PI, geneNames, geneSeqs);
		}
	}
}

void readPhilip(istream &fin){
	int M, L;
	vector<thread> thrds;
	while (fin >> M){
		fin >> L;
		vector<string> geneNames, geneSeqs;
		for (int i = 0; i < M; i++){
			string s, line;
			fin >> s;
			fin >> line;
			geneNames.push_back(s);
			geneSeqs.push_back(line);
		}
		thrds.emplace_back(formatGene, move(geneNames), move(geneSeqs));
		if (thrds.size() == meta.nThreads){
			for (thread &t: thrds) t.join();
			thrds.clear();
		}
	}
	for (thread &t: thrds) t.join();
}

void readFasta(string file){
	ifstream fin(file);
	vector<string> geneNames, geneSeqs;
	string line;
	while (getline(fin, line)){
		if (line[0] == '>'){
			geneNames.push_back(formatName(line));
			geneSeqs.emplace_back();
		}
		else geneSeqs.back() += formatName(line);
	}
	formatGene(geneNames, geneSeqs);
}

string HELP_TEXT = R"V0G0N(-y  take one input in PHYLIP format instead of a list of inputs in FASTA format 
inputList: the path to a file containing a list of paths to input aligned gene files, one file per line
Gene files must be in FASTA format. The header line should be ">Species_Name".
)V0G0N";

int main(int argc, char** argv){
	ARG.setProgramName("asterisk-biallelic", "Accurate Species Tree EstimatoR from Independent Site Kernals");
	ARG.addStringArg('f', "format", "fasta", "Input file type, fasta: a txt file containing a list of FASTA gene files; phylip: a phylip file with all genes");
	ARG.addFlag('F', "formatphylip", "Indicating input file as a PHYLIP file containing all genes (`-f phylip`)", [&](){
		ARG.getStringArg("format") = "phylip";
	}, true);
	string mappingFile;
	meta.initialize(argc, argv, " -y", HELP_TEXT);
	tripInit.numThreads = meta.nThreads;

	bool phylip = ARG.getStringArg("format").compare("fasta");
	
	if (phylip) {
		ifstream fin(ARG.getStringArg("input"));
		readPhilip(fin);
	}
	else {
		ifstream listIn(ARG.getStringArg("input"));
		vector<thread> thrds;
		for (string file; getline(listIn, file);){
			thrds.emplace_back(readFasta, file);
			if (thrds.size() == meta.nThreads){
				for (thread &t: thrds) t.join();
				thrds.clear();
			}
		}
		for (thread &t: thrds) t.join();
	}
	
	cerr << "#Bases: " << tripInit.breakPoints.back() << endl;
	
	score_t score = meta.run().first;
	cerr << "Score: " << (double) score << endl;
	return 0;
}

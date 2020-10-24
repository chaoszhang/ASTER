#include<iostream>
#include<fstream>
#include<unordered_map>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<algorithm>

//#define LARGE_DATA
#ifdef LARGE_DATA
typedef long double score_t;
typedef long long count_t;
#else
typedef double score_t;
typedef int count_t;
#endif

#include "sequence-pair.hpp"
#include "algorithms.hpp"

using namespace std;

TripartitionInitializer tripInit;

int npos = 0;
vector<string> names;
unordered_map<string, int> name2id;
vector<vector<int> > sparse;
int N, n;

string formatName(const string name){
	string res;
	for (char c: name){
		if (c != '>' && c != ' ' && c != '\t') res += c;
	}
	return res;
}

void readFile(istream &fin){
	vector<vector<int> > dist(N, vector<int>(n)), total(N, vector<int>(n));
	vector<vector<bool> > R(N), Y(N);
	vector<pair<score_t, int> > D;
	int id = -1, np = 0, cntR = 0, cntY = 0;
	string line;
	while (getline(fin, line)){
		if (line[0] == '>') {
			id = name2id[formatName(line)];
		}
		else {
			for (char c: line){
				if (c == 'A' || c == 'a' || c == 'G' || c == 'g') { R[id].push_back(1); Y[id].push_back(0); cntR++; }
				if (c == 'C' || c == 'c' || c == 'T' || c == 't') { R[id].push_back(0); Y[id].push_back(1); cntY++; }
				if (c == 'X' || c == 'x' || c == 'N' || c == 'n' || c == '-') { R[id].push_back(0); Y[id].push_back(0); }
			}
			if (R[id].size() > np) np = R[id].size();
		}
	}
	for (id = 0; id < N; id++){
		while (R[id].size() < np) { R[id].push_back(0); Y[id].push_back(0); }
	}
	for (int i = 0; i < N; i++){
		for (int k = 0; k < n; k++){
			int j = sparse[i][k];
			for (int p = 0; p < np; p++){
				dist[i][k] += (R[i][p] && Y[j][p]) || (Y[i][p] && R[j][p]);
				total[i][k] += (R[i][p] || Y[i][p]) && (Y[j][p] || R[j][p]);
			}
		}
	}
	for (int p = 0; p < np; p++){
		score_t P = 0;
		int Q = 0;
		for (int i = 0; i < N; i++){
			for (int k = 0; k < n; k++){
				int j = sparse[i][k];
				if (dist[i][k] > 0 && (R[i][p] || Y[i][p]) && (Y[j][p] || R[j][p])){
					Q++;
					P += ((R[i][p] && Y[j][p]) || (Y[i][p] && R[j][p])) / (score_t) dist[i][k] * total[i][k];
				}
			}
		}
		D.emplace_back(P / Q, p);
	}
	sort(D.begin(), D.end());
	
	for (int r = 0; r + 1 < np; r++){
		int p = D[r].second, q = D[r + 1].second;
		int pcnt[4] = {};
		for (int i = 0; i < N; i++){
			if (R[i][p] && R[i][q]) { tripInit.seq[i].push_back(0); pcnt[0]++; }
			else if (R[i][p] && Y[i][q]) { tripInit.seq[i].push_back(1); pcnt[1]++; }
			else if (Y[i][p] && R[i][q]) { tripInit.seq[i].push_back(2); pcnt[2]++; }
			else if (Y[i][p] && Y[i][q]) { tripInit.seq[i].push_back(3); pcnt[3]++; }
			else tripInit.seq[i].push_back(-1);
		}
		tripInit.pi.push_back(cntR / (score_t) (cntR + cntY));
		tripInit.weight.push_back(pcnt[0] + pcnt[1] + pcnt[2] + pcnt[3] - max({pcnt[0], pcnt[1], pcnt[2], pcnt[3]}) >= 2);
		npos++;
	}
}

void readPhilip(istream &fin){
	string line;
	int nTaxa, L;
	while (fin >> nTaxa){
		vector<vector<int> > dist(N, vector<int>(n)), total(N, vector<int>(n));
		vector<vector<bool> > R(N), Y(N);
		vector<pair<score_t, int> > D;
		int id = -1, np = 0, cntR = 0, cntY = 0;
		fin >> L;
		for (int i = 0; i < nTaxa; i++){
			fin >> line;
			int id = name2id[line];
			getline(fin, line);
			for (char c: line){
				if (c == 'A' || c == 'a' || c == 'G' || c == 'g') { R[id].push_back(1); Y[id].push_back(0); cntR++; }
				if (c == 'C' || c == 'c' || c == 'T' || c == 't') { R[id].push_back(0); Y[id].push_back(1); cntY++; }
				if (c == 'X' || c == 'x' || c == 'N' || c == 'n' || c == '-') { R[id].push_back(0); Y[id].push_back(0); }
			}
			if (R[id].size() > np) np = R[id].size();
		}
		for (id = 0; id < N; id++){
			while (R[id].size() < np) { R[id].push_back(0); Y[id].push_back(0); }
		}
		for (int i = 0; i < N; i++){
			for (int k = 0; k < n; k++){
				int j = sparse[i][k];
				for (int p = 0; p < np; p++){
					dist[i][k] += (R[i][p] && Y[j][p]) || (Y[i][p] && R[j][p]);
					total[i][k] += (R[i][p] || Y[i][p]) && (Y[j][p] || R[j][p]);
				}
			}
		}
		for (int p = 0; p < np; p++){
			score_t P = 0;
			int Q = 0;
			for (int i = 0; i < N; i++){
				for (int k = 0; k < n; k++){
					int j = sparse[i][k];
					if (dist[i][k] > 0 && (R[i][p] || Y[i][p]) && (Y[j][p] || R[j][p])){
						Q++;
						P += ((R[i][p] && Y[j][p]) || (Y[i][p] && R[j][p])) / (score_t) dist[i][k] * total[i][k];
					}
				}
			}
			D.emplace_back(P / Q, p);
		}
		sort(D.begin(), D.end());
		
		for (int r = 0; r + 1 < np; r++){
			int p = D[r].second, q = D[r + 1].second;
			int pcnt[4] = {};
			for (int i = 0; i < N; i++){
				if (R[i][p] && R[i][q]) { tripInit.seq[i].push_back(0); pcnt[0]++; }
				else if (R[i][p] && Y[i][q]) { tripInit.seq[i].push_back(1); pcnt[1]++; }
				else if (Y[i][p] && R[i][q]) { tripInit.seq[i].push_back(2); pcnt[2]++; }
				else if (Y[i][p] && Y[i][q]) { tripInit.seq[i].push_back(3); pcnt[3]++; }
				else tripInit.seq[i].push_back(-1);
			}
			tripInit.pi.push_back(cntR / (score_t) (cntR + cntY));
			tripInit.weight.push_back(pcnt[0] + pcnt[1] + pcnt[2] + pcnt[3] - max({pcnt[0], pcnt[1], pcnt[2], pcnt[3]}) >= 2);
			npos++;
		}
	}
}

string helpText = R"V0G0N(feast [-o oFilePath -r nRound -s nSample -p probability -t nThread -y] inputList
-o  path to output file (default: stdout)
-r  number of total rounds of placements (default: 5)
-s  number of total rounds of subsampling (default: 0)
-p  subsampling probability of keeping each taxon (default: 0.5)
-t  number of threads (default: 1)
-y  take one input in PHYLIP format instead of a list of inputs in FASTA format 
inputList: the path to a file containing a list of paths to input aligned gene files, one file per line
Gene files must be in FASTA format. The header line should be ">Species_Name".
)V0G0N";

int main(int argc, char** argv){
	vector<string> files;
	
	int nThreads = 1, nRounds = 4, nSample = 0;
	bool phylip = false;
	double p = 0.5;
	string outputFile;
	ofstream fileOut;
	if (argc == 1) {cerr << helpText; return 0;}
	for (int i = 1; i < argc; i += 2){
		if (strcmp(argv[i], "-y") == 0) {phylip = true; i--; continue;}
		
		if (strcmp(argv[i], "-o") == 0) outputFile = argv[i + 1];
		if (strcmp(argv[i], "-r") == 0) sscanf(argv[i + 1], "%d", &nRounds);
		if (strcmp(argv[i], "-s") == 0) sscanf(argv[i + 1], "%d", &nSample);
		if (strcmp(argv[i], "-p") == 0) sscanf(argv[i + 1], "%lf", &p);
		if (strcmp(argv[i], "-t") == 0) sscanf(argv[i + 1], "%d", &nThreads);
		if (strcmp(argv[i], "-h") == 0) {cerr << helpText; return 0;}
	}
	ostream &fout = (outputFile == "") ? cout : fileOut;
	if (outputFile != "") fileOut.open(outputFile);
	
	if (nRounds < nThreads){
		tripInit.nThreads = nThreads / nRounds;
		nThreads = nRounds;
	}
	
	if (phylip) {
		{
			ifstream fin(argv[argc - 1]);
			int n, L;
			string line;
			while (fin >> n){
				fin >> L;
				for (int i = 0; i < n; i++){
					string s;
					fin >> s;
					getline(fin, line);
					if(name2id.count(s) == 0) {
						name2id[s] = names.size();
						names.push_back(s);
						tripInit.seq.push_back({});
					}
				}
			}
		}
		N = names.size();
		n = sqrt(N);
		for (int i = 0; i < N; i++){
			sparse.emplace_back();
			for (int j = 0; j < 5 * n; j++){
				sparse[i].push_back(rand() % N);
			}
		}
		{
			ifstream fin(argv[argc - 1]);
			readPhilip(fin);
		}
	}
	else {
		ifstream listIn(argv[argc - 1]);
		for (string file; getline(listIn, file);){
			files.push_back(file);
		}
		
		for (string file: files){
			ifstream fin(file);
			string line;
			while (getline(fin, line)){
				if (line[0] != '>') continue;
				string s = formatName(line);
				if(name2id.count(s) == 0) {
					name2id[s] = names.size();
					names.push_back(s);
					tripInit.seq.push_back({});
				}
			}
		}
		N = names.size();
		n = sqrt(N);
		for (int i = 0; i < N; i++){
			sparse.emplace_back();
			for (int j = 0; j < 5 * n; j++){
				sparse[i].push_back(rand() % N);
			}
		}
		for (string file: files){
			ifstream fin(file);
			readFile(fin);
		}
	}
	
	for (int i = 0; i < names.size(); i++){
		tripInit.weightHelper.push_back(0);
	}
	tripInit.weightHelper.push_back(0);
	
	tripInit.npos = npos;
	cerr << "#Species: " << names.size() << endl;
	cerr << "#Bases: " << npos << endl;
	cerr << "#Rounds: " << nRounds << endl;
	cerr << "#Samples: " << nSample << endl;
	cerr << "#Threads: " << nThreads << "x" << tripInit.nThreads << endl;
	cerr << "p = " << p << endl;
	
	ConstrainedOptimizationAlgorithm alg(names.size(), tripInit, names);
	auto res = alg.run(nRounds, nThreads);
	cerr << "Score: " << (double) res.first << endl;
	cerr << res.second << endl;
	
	res = alg.run(nSample, nThreads, p);
	cerr << "Score: " << (double) res.first << endl;
	fout << res.second << endl;
	return 0;
}

#include<iostream>
#include<fstream>
#include<unordered_map>
#include<cstdio>
#include<cstdlib>
#include<cstring>

//#define LARGE_DATA
#ifdef LARGE_DATA
typedef long double score_t;
typedef long long count_t;
#else
typedef double score_t;
typedef int count_t;
#endif

#include "sequence.hpp"
#include "algorithms.hpp"

using namespace std;

int npos = 0;
vector<array<score_t, 4> > pi;
vector<vector<int> > seq;
vector<string> names;
unordered_map<string, int> name2id;

string formatName(const string name){
	string res;
	for (char c: name){
		if (c != '>' && c != ' ' && c != '\t') res += c;
	}
	return res;
}

void readFile(istream &fin){
	int cnt[4] = {}, id = -1;
	string line;
	while (getline(fin, line)){
		if (line[0] == '>') {
			if (id != -1 && npos < seq[id].size()) npos = seq[id].size();
			id = name2id[formatName(line)];
		}
		else{
			for (char c: line){
				if (c == 'A' || c == 'a') { seq[id].push_back(0); cnt[0]++; }
				if (c == 'C' || c == 'c') { seq[id].push_back(1); cnt[1]++; }
				if (c == 'G' || c == 'g') { seq[id].push_back(2); cnt[2]++; }
				if (c == 'T' || c == 't') { seq[id].push_back(3); cnt[3]++; }
				if (c == 'X' || c == 'x' || c == 'N' || c == 'n' || c == '-') { seq[id].push_back(-1); }
			}
		}
	}
	if (npos < seq[id].size()) npos = seq[id].size();
	for (id = 0; id < seq.size(); id++){
		while (seq[id].size() < npos) seq[id].push_back(-1);
	}
	score_t cntsum = cnt[0] + cnt[1] + cnt[2] + cnt[3];
	while (pi.size() < npos) pi.push_back({cnt[0] / cntsum, cnt[1] / cntsum, cnt[2] / cntsum, cnt[3] / cntsum});
}

void readPhilip(istream &fin){
	int cnt[4] = {}, id = -1;
	string line;
	int nTaxa, L;
	while (fin >> nTaxa){
		fin >> L;
		npos += L;
		for (int i = 0; i < nTaxa; i++){
			fin >> line;
			int id = name2id[line];
			getline(fin, line);
			for (char c: line){
				if (c == 'A' || c == 'a') { seq[id].push_back(0); cnt[0]++; }
				if (c == 'C' || c == 'c') { seq[id].push_back(1); cnt[1]++; }
				if (c == 'G' || c == 'g') { seq[id].push_back(2); cnt[2]++; }
				if (c == 'T' || c == 't') { seq[id].push_back(3); cnt[3]++; }
				if (c == 'X' || c == 'x' || c == 'N' || c == 'n' || c == '-') { seq[id].push_back(-1); }
			}
			for (id = 0; id < seq.size(); id++){
				while (seq[id].size() < npos) seq[id].push_back(-1);
			}
			score_t cntsum = cnt[0] + cnt[1] + cnt[2] + cnt[3];
			while (pi.size() < npos) pi.push_back({cnt[0] / cntsum, cnt[1] / cntsum, cnt[2] / cntsum, cnt[3] / cntsum});
		}
	}
}

string helpText = R"V0G0N(faster [-o oFilePath -r nRound -s nSample -p probability -t nThread -y] inputList
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
	
	int nThreads = 1, nRounds = 5, nSample = 0;
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
						seq.push_back({});
					}
				}
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
					seq.push_back({});
				}
			}
		}
		for (string file: files){
			ifstream fin(file);
			readFile(fin);
		}
	}
	
	
	TripartitionInitializer tripInit;
	tripInit.npos = npos;
	tripInit.pi = pi;
	tripInit.seq = seq;
	cerr << "#Species: " << names.size() << endl;
	cerr << "#Bases: " << npos << endl;
	cerr << "#Rounds: " << nRounds << endl;
	cerr << "#Samples: " << nSample << endl;
	cerr << "#Threads: " << nThreads << endl;
	cerr << "p = " << p << endl;
	
	ConstrainedOptimizationAlgorithm alg(names.size(), tripInit, names);
	auto res = alg.run(nRounds, nThreads);
	cerr << "Score: " << res.first << endl;
	cerr << res.second << endl;
	
	res = alg.run(nRounds, nThreads, 1);
	cerr << "Score: " << res.first << endl;
	fout << res.second << endl;
	return 0;
}

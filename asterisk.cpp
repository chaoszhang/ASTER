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

MetaAlgorithm meta;
TripartitionInitializer &tripInit = meta.tripInit;

int &npos = meta.tripInit.npos;
vector<string> &names = meta.names;
unordered_map<string, int> &name2id = meta.name2id;
vector<vector<int> > sparse;
int N, n;

string formatName(const string name){
	string res;
	for (char c: name){
		if (c != '>' && c != ' ' && c != '\t') res += c;
	}
	return res;
}

void readFile(istream &fin, const char CR1, const char CR2, const char CY1, const char CY2, const char CY3, const char cr1, const char cr2, const char cy1, const char cy2, const char cy3){
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
				if (c == CR1 || c == cr1 || c == CR2 || c == cr2) { R[id].push_back(1); Y[id].push_back(0); cntR++; }
				if (c == CY1 || c == cy1 || c == CY2 || c == cy2 || c == CY3 || c == cy3) { R[id].push_back(0); Y[id].push_back(1); cntY++; }
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

void readPhilip(istream &fin, const char CR1, const char CR2, const char CY1, const char CY2, const char CY3, const char cr1, const char cr2, const char cy1, const char cy2, const char cy3){
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
				if (c == CR1 || c == cr1 || c == CR2 || c == cr2) { R[id].push_back(1); Y[id].push_back(0); cntR++; }
				if (c == CY1 || c == cy1 || c == CY2 || c == cy2 || c == CY3 || c == cy3) { R[id].push_back(0); Y[id].push_back(1); cntY++; }
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

string HELP_TEXT = R"V0G0N(-y  take one input in PHYLIP format instead of a list of inputs in FASTA format 
inputList: the path to a file containing a list of paths to input aligned gene files, one file per line
Gene files must be in FASTA format. The header line should be ">Species_Name".
)V0G0N";

int main(int argc, char** argv){
	vector<string> files;
	bool phylip = false;
	string mappingFile;
	meta.initialize(argc, argv, " -y", HELP_TEXT);
	
	for (int i = 1; i < argc; i += 2){
		if (strcmp(argv[i], "-y") == 0) {phylip = true; i--; continue;}
	}
	
	if (phylip) {
		{
			ifstream fin(argv[argc - 1]);
			int M, L;
			string line;
			while (fin >> M){
				fin >> L;
				for (int i = 0; i < M; i++){
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
		n = 2 * log(N) * sqrt(N);
		for (int i = 0; i < N; i++){
			sparse.emplace_back();
			for (int j = 0; j < n; j++){
				sparse[i].push_back(rand() % N);
			}
		}
		{ ifstream fin(argv[argc - 1]); readPhilip(fin, 'A', 'G', 'C', 'T', 0, 'a', 'g', 'c', 't', 0); }
		{ ifstream fin(argv[argc - 1]); readPhilip(fin, 'A', 'C', 'G', 'T', 0, 'a', 'c', 'g', 't', 0); }
		{ ifstream fin(argv[argc - 1]); readPhilip(fin, 'A', 'T', 'C', 'G', 0, 'a', 't', 'c', 'g', 0); }
		{ ifstream fin(argv[argc - 1]); readPhilip(fin, 'A', 0, 'T', 'C', 'G', 'a', 0, 't', 'c', 'g'); }
		{ ifstream fin(argv[argc - 1]); readPhilip(fin, 'T', 0, 'A', 'C', 'G', 't', 0, 'a', 'c', 'g'); }
		{ ifstream fin(argv[argc - 1]); readPhilip(fin, 'C', 0, 'A', 'T', 'G', 'c', 0, 'a', 't', 'g'); }
		{ ifstream fin(argv[argc - 1]); readPhilip(fin, 'G', 0, 'A', 'T', 'C', 'g', 0, 'a', 't', 'c'); }
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
		n = 2 * log(N) * sqrt(N);
		for (int i = 0; i < N; i++){
			sparse.emplace_back();
			for (int j = 0; j < n; j++){
				sparse[i].push_back(rand() % N);
			}
		}
		for (string file: files){
			{ ifstream fin(file); readFile(fin, 'A', 'G', 'C', 'T', 0, 'a', 'g', 'c', 't', 0); }
			{ ifstream fin(file); readFile(fin, 'A', 'C', 'G', 'T', 0, 'a', 'c', 'g', 't', 0); }
			{ ifstream fin(file); readFile(fin, 'A', 'T', 'C', 'G', 0, 'a', 't', 'c', 'g', 0); }
			{ ifstream fin(file); readFile(fin, 'A', 0, 'T', 'C', 'G', 'a', 0, 't', 'c', 'g'); }
			{ ifstream fin(file); readFile(fin, 'T', 0, 'A', 'C', 'G', 't', 0, 'a', 'c', 'g'); }
			{ ifstream fin(file); readFile(fin, 'C', 0, 'A', 'T', 'G', 'c', 0, 'a', 't', 'g'); }
			{ ifstream fin(file); readFile(fin, 'G', 0, 'A', 'T', 'C', 'g', 0, 'a', 't', 'c'); }
		}
	}
	
	cerr << "#Bases: " << npos << endl;
	
	score_t score = meta.run().first;
	cerr << "Score: " << (double) score << endl;
	return 0;
}

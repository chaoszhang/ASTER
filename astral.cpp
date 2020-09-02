#include<iostream>
#include<fstream>
#include<unordered_map>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<string>

using namespace std;

//#define LARGE_DATA
#ifdef LARGE_DATA
typedef __int128 score_t;

string to_string(const __int128 x){
	return to_string((double) x);
}

ostream& operator<<(ostream& cout, __int128 x){
	return cout << to_string(x);
}
#else
typedef long long score_t;
#endif

#include "genetree.hpp"
#include "algorithms.hpp"


TripartitionInitializer tripInit;

unordered_map<string, string> leafname_mapping;
string TEXT;
int nodecnt = 0;
int pos = 0;
vector<string> names;
unordered_map<string, int> name2id;


string MAPPING(int begin, int end){
	string s;
	for (int i = begin; i < end && TEXT[i] != ':'; i++){
		if (TEXT[i] != '\"' && TEXT[i] != '\'') s += TEXT[i];
	}
	if (leafname_mapping.count(s)) return leafname_mapping[s];
	else return s;
}

void parseTaxa(){
	int i = pos;
	while (TEXT[pos] != '(' && TEXT[pos] != ',') pos++;
	if (TEXT[pos] == '(') {
		pos++;
		parseTaxa();
		while (TEXT[pos] != ',') pos++;
	}
	else {
		string name = MAPPING(i, pos);
		if (name2id.count(name) == 0){
			name2id[name] = names.size();
			names.push_back(name);
		}
	}
	while (TEXT[pos] != ')'){
		i = ++pos;
		while (TEXT[pos] != ')' && TEXT[pos] != '(' && TEXT[pos] != ',') pos++;
		if (TEXT[pos] == '(') {
			pos++;
			parseTaxa();
			while (TEXT[pos] != ',' && TEXT[pos] != ')') pos++;
		}
		else {
			string name = MAPPING(i, pos);
			if (name2id.count(name) == 0){
				name2id[name] = names.size();
				names.push_back(name);
			}
		}
	}
	pos++;
}

int parse(vector<pair<int, int> > &parentChild){
	int i = pos;
	int cur = nodecnt++;
	while (TEXT[pos] != '(' && TEXT[pos] != ',') pos++;
	if (TEXT[pos] == '(') {
		pos++;
		parentChild.emplace_back(cur, parse(parentChild));
		while (TEXT[pos] != ',') pos++;
	}
	else parentChild.emplace_back(cur, name2id[MAPPING(i, pos)]);
	while (TEXT[pos] != ')'){
		i = ++pos;
		while (TEXT[pos] != ')' && TEXT[pos] != '(' && TEXT[pos] != ',') pos++;
		if (TEXT[pos] == '(') {
			pos++;
			parentChild.emplace_back(cur, parse(parentChild));
			while (TEXT[pos] != ',' && TEXT[pos] != ')') pos++;
		}
		else parentChild.emplace_back(cur, name2id[MAPPING(i, pos)]);
	}
	pos++;
	return cur;
}

void readInputTrees(string input, string mapping) {
	if (mapping != ""){
		ifstream fmap(mapping);
		string gname, sname;
		while (fmap >> gname){
			fmap >> sname;
			leafname_mapping[gname] = sname;
		}
	}
	ifstream fin(input);
	string line;
	while (getline(fin, line)) TEXT += line;
	
	while (pos < TEXT.size()){
		while (pos < TEXT.size() && TEXT[pos] != '(') pos++;
		if (pos < TEXT.size()) {
			pos++;
			parseTaxa();
		}
	}
	tripInit.nTaxa = names.size();
	pos = 0;
	while (pos < TEXT.size()){
		while (pos < TEXT.size() && TEXT[pos] != '(') pos++;
		if (pos < TEXT.size()) {
			pos++;
			tripInit.parentChild.emplace_back();
			nodecnt = names.size();
			tripInit.roots.push_back(parse(tripInit.parentChild.back()));
		}
	}
}

string helpText = R"V0G0N(astral_feast [-o oFilePath -r nRound -s nSample -p probability -t nThread -a taxonNameMaps] inputGeneTrees
-o  path to output file (default: stdout)
-r  number of total rounds of placements (default: 5)
-s  number of total rounds of subsampling (default: 0)
-p  subsampling probability of keeping each taxon (default: 0.5)
-t  number of threads (default: 1)
-a  a list of gene name to taxon name maps, each line contains one gene name followed by one taxon name separated by a space or tab 
inputGeneTrees: the path to a file containing all gene trees in Newick format
)V0G0N";

int main(int argc, char** argv){
	int nThreads = 1, nRounds = 4, nSample = 0;
	bool phylip = false;
	double p = 0.5;
	string outputFile, mappingFile;
	ofstream fileOut;
	if (argc == 1) {cerr << helpText; return 0;}
	for (int i = 1; i < argc; i += 2){
		if (strcmp(argv[i], "-a") == 0) mappingFile = argv[i + 1];
		
		if (strcmp(argv[i], "-o") == 0) outputFile = argv[i + 1];
		if (strcmp(argv[i], "-r") == 0) sscanf(argv[i + 1], "%d", &nRounds);
		if (strcmp(argv[i], "-s") == 0) sscanf(argv[i + 1], "%d", &nSample);
		if (strcmp(argv[i], "-p") == 0) sscanf(argv[i + 1], "%lf", &p);
		if (strcmp(argv[i], "-t") == 0) sscanf(argv[i + 1], "%d", &nThreads);
		if (strcmp(argv[i], "-h") == 0) {cerr << helpText; return 0;}
	}
	ostream &fout = (outputFile == "") ? cout : fileOut;
	if (outputFile != "") fileOut.open(outputFile);
	readInputTrees(argv[argc - 1], mappingFile);
	
	/*
	for (auto e: tripInit.parentChild[0]){
		cerr << e.first << " ";
		if (e.second < tripInit.nTaxa) cerr << names[e.second] << endl; else cerr << e.second << endl;
	}
	*/
	
	
	cerr << "#Species: " << names.size() << endl;
	cerr << "#Genetrees: " << tripInit.parentChild.size() << endl;
	cerr << "#Rounds: " << nRounds << endl;
	cerr << "#Samples: " << nSample << endl;
	cerr << "#Threads: " << nThreads << endl;
	cerr << "p = " << p << endl;
	
	ConstrainedOptimizationAlgorithm alg(names.size(), tripInit, names);
	auto res = alg.run(nRounds, nThreads);
	cerr << "Score: " << res.first/2 << endl;
	cerr << res.second << endl;
	
	res = alg.run(nSample, nThreads, p);
	cerr << "Score: " << res.first/2 << endl;
	fout << res.second << endl;
	
	//cerr << alg.printOptimalTreeWithScore() << endl;
	return 0;
}

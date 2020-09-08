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
typedef long double score_t;

string to_string(long double x){
	return to_string((double) x);
}

score_t from_string(const string s){
	return stold(s);
}

ostream& operator<<(ostream& cout, long double x){
	return cout << to_string(x);
}
#else
typedef double score_t;

score_t from_string(const string s){
	return stod(s);
}
#endif

#include "genetreewithlength.hpp"
#include "algorithms.hpp"


TripartitionInitializer tripInit;

unordered_map<string, string> leafname_mapping;
string TEXT;
int pos = 0;
int K = 0;
vector<string> names;
unordered_map<string, int> name2id;


int MAPPING(int begin, int end){
	string s;
	for (int i = begin; i < end && TEXT[i] != ':'; i++){
		if (TEXT[i] != '\"' && TEXT[i] != '\'') s += TEXT[i];
	}
	if (leafname_mapping.count(s)) s = leafname_mapping[s];
	if (name2id.count(s) == 0){
		name2id[s] = names.size();
		names.push_back(s);
		tripInit.leafParent.emplace_back();
	}
	return name2id[s];
}

score_t WEIGHT(int begin, int end){
	int i = begin;
	while (i < end && TEXT[i] != ':') i++;
	if (i == end) return 1;
	else return exp(-from_string(TEXT.substr(i + 1, end - i - 1)));
}

void parse(int parent = -1, bool isLeft = true){
	int cur = tripInit.nodes.size();
	tripInit.nodes.emplace_back();
	tripInit.nodes[cur].up = parent;
	if (parent != -1 && isLeft) tripInit.nodes[parent].small = cur;
	if (parent != -1 && !isLeft) tripInit.nodes[parent].large = cur;
	
	if (TEXT[pos] == '(') { //cerr << TEXT[pos];
		pos++;
		parse(cur, true); //cerr << TEXT[pos];
		pos++;
		parse(cur, false); //cerr << TEXT[pos];
		while (TEXT[pos] != ')'){
			int next = tripInit.nodes.size();
			tripInit.nodes.emplace_back();
			tripInit.nodes[cur].up = next;
			tripInit.nodes[next].small = cur;
			cur = next;
			tripInit.nodes[cur].up = parent;
			if (parent != -1 && isLeft) tripInit.nodes[parent].small = cur;
			if (parent != -1 && !isLeft) tripInit.nodes[parent].large = cur;
			pos++;
			parse(cur, false); //cerr << TEXT[pos];
		}
		int i = ++pos;
		while (TEXT[pos] != ')' && TEXT[pos] != ',' && TEXT[pos] != ';') pos++;
		tripInit.nodes[cur].weight = WEIGHT(i, pos); //cerr << ":" << tripInit.nodes[cur].weight;
	} 
	else {
		int i = pos;
		while (TEXT[pos] != ')' && TEXT[pos] != ',') pos++;
		tripInit.leafParent[MAPPING(i, pos)].push_back(cur);
		tripInit.nodes[cur].weight = WEIGHT(i, pos); //cerr << "'"<< names[MAPPING(i, pos)] << "':" << tripInit.nodes[cur].weight;
	}
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
			parse(); //cerr << endl;
			K++;
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
	cerr << "#Genetrees: " << K << endl;
	cerr << "#Rounds: " << nRounds << endl;
	cerr << "#Samples: " << nSample << endl;
	cerr << "#Threads: " << nThreads << endl;
	cerr << "p = " << p << endl;
	
	ConstrainedOptimizationAlgorithm alg(names.size(), tripInit, names);
	auto res = alg.run(nRounds, nThreads);
	cerr << "Score: " << res.first << endl;
	cerr << res.second << endl;
	
	res = alg.run(nSample, nThreads, p);
	cerr << "Score: " << res.first << endl;
	fout << res.second << endl;
	
	//cerr << alg.printOptimalTreeWithScore() << endl;
	return 0;
}

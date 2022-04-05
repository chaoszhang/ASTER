#define DRIVER_VERSION "1"

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

#include "argparser.hpp"
#include "genetreewithlength.hpp"
#include "algorithms.hpp"


MetaAlgorithm meta;
TripartitionInitializer &tripInit = meta.tripInit;
vector<TripartitionInitializer> &batchInit = meta.batchInit;

unordered_map<string, unordered_set<string> > reverse_mapping;
unordered_map<string, string> leafname_mapping;
string TEXT;
int pos = 0;
int K = 0;
int part = 0, iBatch = 0;
vector<string> &names = meta.names;
unordered_map<string, int> &name2id = meta.name2id;
vector<int> &nameCnts = meta.tripInit.nameCnts;
double lengthFactor = -1;

int MAPPING(int begin, int end){
	string s;
	for (int i = begin; i < end && TEXT[i] != ':'; i++){
		if (TEXT[i] != '\"' && TEXT[i] != '\'') s += TEXT[i];
	}
	if (leafname_mapping.count(s)) {
		reverse_mapping[leafname_mapping[s]].insert(s);
		s = leafname_mapping[s];
	}
	else{
		reverse_mapping[s].insert(s);
	}
	if (name2id.count(s) == 0){
		name2id[s] = names.size();
		names.push_back(s);
		for (int i = 0; i < meta.nThread2; i++){
			tripInit.leafParent[i].emplace_back();
		}
	}
	return name2id[s];
}

score_t WEIGHT(int begin, int end){
	int i = begin;
	while (i < end && TEXT[i] != ':') i++;
	if (i == end) return 1;
	else return exp(lengthFactor * from_string(TEXT.substr(i + 1, end - i - 1)));
}

void parse(int parent = -1, bool isLeft = true){
	int cur = tripInit.nodes[part].size();
	tripInit.nodes[part].emplace_back();
	tripInit.nodes[part][cur].up = parent;
	if (parent != -1 && isLeft) tripInit.nodes[part][parent].small = cur;
	if (parent != -1 && !isLeft) tripInit.nodes[part][parent].large = cur;
	
	if (TEXT[pos] == '(') { 
		pos++;
		parse(cur, true); 
		pos++;
		parse(cur, false); 
		vector<int> lst;
		lst.push_back(cur);
		tripInit.nodes[part][cur].weight = 1;
		while (TEXT[pos] != ')'){
			int left = lst[rand() % lst.size()];
			int up = tripInit.nodes[part].size();
			tripInit.nodes[part].emplace_back();
			lst.push_back(up);
			if (cur == left) cur = up;
			tripInit.nodes[part][up].weight = 1;
			int g = tripInit.nodes[part][left].up;
			if (g != -1){
				if (tripInit.nodes[part][g].small == left) tripInit.nodes[part][g].small = up;
				else tripInit.nodes[part][g].large = up;
			}
			tripInit.nodes[part][up].up = g;
			tripInit.nodes[part][left].up = up; 
			tripInit.nodes[part][up].small = left;
			pos++;
			parse(up, false);
		}
		int i = ++pos;
		while (TEXT[pos] != ')' && TEXT[pos] != ',' && TEXT[pos] != ';') pos++;
		tripInit.nodes[part][cur].weight = WEIGHT(i, pos);
	} 
	else {
		int i = pos;
		while (TEXT[pos] != ')' && TEXT[pos] != ',') pos++;
		tripInit.leafParent[part][MAPPING(i, pos)].push_back(cur);
		tripInit.nodes[part][cur].weight = WEIGHT(i, pos);
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
			part = K % tripInit.nodes.size();
			parse();
			K++;
		}
	}
}

string HELP = " -a taxonNameMaps";
string HELP_TEXT = R"V0G0N(-a  a list of gene name to taxon name maps, each line contains one gene name followed by one taxon name separated by a space or tab 
inputGeneTrees: the path to a file containing all gene trees in Newick format
)V0G0N";

int main(int argc, char** argv){
	ARG.setProgramName("astral-lengthweighted", "Weighted ASTRAL by Branch Length");
	ARG.addStringArg('a', "mapping", "", "A list of gene name to taxon name maps, each line contains one gene name followed by one taxon name separated by a space or tab");
	ARG.addDoubleArg('w', "weight", -1, "Weight factor of total terminal branch lengths");
	/*ARG.addFlag('H', "hgt", "Use this preset when you are sure that horizontal gene transfer is the main source of discordance. (`-w 1`)", [&](){
		ARG.getDoubleArg("weight") = 1;
	}, true);*/
	int dupType = 1;
	string mappingFile;
	meta.initialize(argc, argv, HELP, HELP_TEXT);
	mappingFile = ARG.getStringArg("mapping");
	lengthFactor = ARG.getDoubleArg("weight");

	for (int i = 0; i < meta.nThread2; i++){
		tripInit.nodes.emplace_back();
		tripInit.leafParent.emplace_back();
	}
	for (int i = 0; i < meta.nBatch; i++){
		batchInit[i].nodes.emplace_back();
		batchInit[i].leafParent.emplace_back();
	}
	readInputTrees(ARG.getStringArg("input"), mappingFile);
	
	if (dupType == 2){
		for (string s: names){
			nameCnts.push_back(reverse_mapping[s].size());
		}
	}
	else {
		for (string s: names){
			nameCnts.push_back(1);
		}
	}
	cerr << "#Genetrees: " << K << endl;
	
	score_t score = meta.run().first;
	//cerr << "Score: " << (double) score << endl;
	fprintf(stderr, "Score: %.10lg\n", (double) score);
	return 0;
}

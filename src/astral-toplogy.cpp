#define DRIVER_VERSION "0"

/* CHANGE LOG
 */

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

#include "argparser.hpp"
#include "genetree.hpp"
#include "algorithms.hpp"

MetaAlgorithm meta;
TripartitionInitializer &tripInit = meta.tripInit;
vector<TripartitionInitializer> &batchInit = meta.batchInit;

unordered_map<string, string> leafname_mapping;
string TEXT;
int nodecnt = 0;
int pos = 0;
int K = 0;
int part = 0, iBatch = 0;
vector<string> &names = meta.names;
unordered_map<string, int> &name2id = meta.name2id;


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
	while (TEXT[pos] != '(' && TEXT[pos] != ',' && TEXT[pos] != ')') pos++;
	
	if (TEXT[pos] == '(') {
		while (TEXT[pos] != ')'){
			pos++;
			parseTaxa();
			while (TEXT[pos] != ')' && TEXT[pos] != ',') pos++;
		}
		pos++;
	}
	else {
		string name = MAPPING(i, pos);
		if (name2id.count(name) == 0){
			name2id[name] = names.size();
			names.push_back(name);
		}
	}
}

void parse(int parent = -1){
	int i = pos, cur;
	while (TEXT[pos] != '(' && TEXT[pos] != ',' && TEXT[pos] != ')') pos++;
	
	if (TEXT[pos] == '(') {
		cur = nodecnt++;
		while (TEXT[pos] != ')'){
			pos++;
			parse(cur);
			while (TEXT[pos] != ')' && TEXT[pos] != ',') pos++;
		}
		pos++;
	}
	else cur = name2id[MAPPING(i, pos)];
	if (parent != -1) tripInit.parentChild[part].back().push_back({parent, cur});
	else tripInit.roots[part].push_back(cur);
	if (parent != -1) batchInit[iBatch].parentChild[0].back().push_back({parent, cur});
	else batchInit[iBatch].roots[0].push_back(cur);
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
		parseTaxa();
		while (pos < TEXT.size() && TEXT[pos] != '(') pos++;
	}
	tripInit.nTaxa = names.size();
	for (int i = 0; i < meta.nBatch; i++){
		batchInit[i].nTaxa = names.size();
	}
	pos = 0;
	while (pos < TEXT.size()){
		part = K % tripInit.roots.size();
		tripInit.parentChild[part].emplace_back();
		iBatch = K % meta.nBatch;
		batchInit[iBatch].parentChild[0].emplace_back();
		nodecnt = names.size();
		parse();
		while (pos < TEXT.size() && TEXT[pos] != '(') pos++;
		K++;
	}
}

string HELP_TEXT = R"V0G0N(-a  a list of gene name to taxon name maps, each line contains one gene name followed by one taxon name separated by a space or tab 
inputGeneTrees: the path to a file containing all gene trees in Newick format
)V0G0N";

int main(int argc, char** argv){
	ARG.setProgramName("astral", "Accurate Species TRee ALgorithm");
	ARG.addStringArg('a', "mapping", "", "A list of gene name to taxon name maps, each line contains one gene name followed by one taxon name separated by a space or tab");
	
	string mappingFile;
	meta.initialize(argc, argv, " -a taxonNameMaps", HELP_TEXT);
	mappingFile = ARG.getStringArg("mapping");
	
	for (int i = 0; i < meta.nThread2; i++){
		tripInit.roots.emplace_back();
		tripInit.parentChild.emplace_back();
	}
	for (int i = 0; i < meta.nBatch; i++){
		batchInit[i].roots.emplace_back();
		batchInit[i].parentChild.emplace_back();
	}
	readInputTrees(ARG.getStringArg("input"), mappingFile);
	
	cerr << "#Genetrees: " << K << endl;
	
	score_t score = meta.run().first / 2;
	cerr << "Score: " << score << endl;
	return 0;
}

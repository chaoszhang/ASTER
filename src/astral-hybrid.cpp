#define DRIVER_VERSION "7"

/* CHANGE LOG
 * 7: bug fix
 * 6: auto-detect support type
 * 5: merge all versions of weighted astral
 * 4: add weighting by tree and other options
 */

#define ROOTING

#include<iostream>
#include<fstream>
#include<unordered_map>
#include<unordered_set>
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
#include "genetreewithhybridweight.hpp"
#include "algorithms.hpp"


MetaAlgorithm meta;
TripartitionInitializer &tripInit = meta.tripInit;

unordered_map<string, unordered_set<string> > reverse_mapping;
unordered_map<string, string> leafname_mapping;
string TEXT;
int pos = 0;
int K = 0;
int part = 0, iBatch = 0;
vector<string> &names = meta.names;
unordered_map<string, int> &name2id = meta.name2id;
vector<int> &nameCnts = meta.tripInit.nameCnts;
score_t maxv = 100, minv = 0, defaultv = 0;
double lengthFactor = -1;
bool useSupport = true, useLength = true;

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
		for (int i = 0; i < meta.nThreads; i++){
			tripInit.leafParent[i].emplace_back();
		}
	}
	return name2id[s];
}

score_t WEIGHT_S(int begin, int end){
	if (!useSupport) return 1;
	int i = begin;
	while (i < end && TEXT[i] != ':') i++;
	if (i == begin) return max((score_t)0.0, (defaultv - minv) / (maxv - minv));
	else return max((score_t)0.0, (from_string(TEXT.substr(begin, i - begin)) - minv) / (maxv - minv));
}

score_t DEFAULT_WEIGHT_S(){
	if (!useSupport) return 0;
	else return max((score_t)0.0, (defaultv - minv) / (maxv - minv));
}
score_t WEIGHT_L(int begin, int end){
	if (!useLength) return 1;
	int i = begin;
	while (i < end && TEXT[i] != ':') i++;
	if (i == end) return 1;
	else return exp(lengthFactor * from_string(TEXT.substr(i + 1, end - i - 1)));
}

void parse(score_t treeweight, int parent = -1, bool isLeft = true){
	int cur = tripInit.nodes[part].size();
	tripInit.nodes[part].emplace_back();
	tripInit.nodes[part][cur].up = parent;
	if (parent != -1 && isLeft) tripInit.nodes[part][parent].small = cur;
	if (parent != -1 && !isLeft) tripInit.nodes[part][parent].large = cur;
	
	if (TEXT[pos] == '(') { 
		pos++;
		parse(treeweight, cur, true); 
		pos++;
		parse(treeweight, cur, false);
		vector<int> lst;
		lst.push_back(cur);
		tripInit.nodes[part][cur].weight = DEFAULT_WEIGHT_S();
		tripInit.nodes[part][cur].length = 1;
		tripInit.nodes[part][cur].treeweight = treeweight;
		while (TEXT[pos] != ')'){
			int left = lst[rand() % lst.size()];
			int up = tripInit.nodes[part].size();
			tripInit.nodes[part].emplace_back();
			lst.push_back(up);
			if (cur == left) cur = up;
			tripInit.nodes[part][up].weight = DEFAULT_WEIGHT_S();
			tripInit.nodes[part][up].length = 1;
			tripInit.nodes[part][up].treeweight = treeweight;
			int g = tripInit.nodes[part][left].up;
			if (g != -1){
				if (tripInit.nodes[part][g].small == left) tripInit.nodes[part][g].small = up;
				else tripInit.nodes[part][g].large = up;
			}
			tripInit.nodes[part][up].up = g;
			tripInit.nodes[part][left].up = up; 
			tripInit.nodes[part][up].small = left;
			pos++;
			parse(treeweight, up, false);
		}
		int i = ++pos;
		while (TEXT[pos] != ')' && TEXT[pos] != ',' && TEXT[pos] != ';') pos++;
		tripInit.nodes[part][cur].weight = WEIGHT_S(i, pos);
		tripInit.nodes[part][cur].length = WEIGHT_L(i, pos);
		tripInit.nodes[part][cur].treeweight = treeweight;
	} 
	else {
		int i = pos;
		while (TEXT[pos] != ')' && TEXT[pos] != ',') pos++;
		tripInit.leafParent[part][MAPPING(i, pos)].push_back(cur);
		tripInit.nodes[part][cur].weight = 1;
		tripInit.nodes[part][cur].length = WEIGHT_L(i, pos);
		tripInit.nodes[part][cur].treeweight = treeweight;
	}
}

void readInputTrees(string input, string mapping, string treeweights) {
	if (mapping != ""){
		ifstream fmap(mapping);
		string gname, sname;
		while (fmap >> gname){
			fmap >> sname;
			leafname_mapping[gname] = sname;
		}
	}
	ifstream fin(input), ftw;
	if (treeweights != "") ftw.open(treeweights);

	string line;
	while (getline(fin, line)) TEXT += line;
	while (pos < TEXT.size()){
		while (pos < TEXT.size() && TEXT[pos] != '(') pos++;
		if (pos < TEXT.size()) {
			part = K % tripInit.nodes.size();
			score_t tw = 1;
			if (treeweights != "") ftw >> tw;
			parse(tw);
			K++;
		}
	}
}

void DETECT(int begin, int end){
	int i = begin;
	while (i < end && TEXT[i] != ':') i++;
	if (i == begin) return;
	score_t s = from_string(TEXT.substr(begin, i - begin));
	if (s > maxv) maxv = s;
	if (s < minv) minv = s;
}

void detectSupport(){
	if (TEXT[pos] == '(') { 
		pos++;
		detectSupport(); 
		pos++;
		detectSupport();
		while (TEXT[pos] != ')') {
			pos++;
			detectSupport();
		}
		int i = ++pos;
		while (TEXT[pos] != ')' && TEXT[pos] != ',' && TEXT[pos] != ';') pos++;
		DETECT(i, pos);
	} 
	else while (TEXT[pos] != ')' && TEXT[pos] != ',') pos++;
}

void detectSupportType(string input) {
	LOG << "Auto-detecting support value type ...\n";
	ifstream fin(input);
	string line;
	maxv = 0;
	minv = 100;

	while (getline(fin, line)) TEXT += line;
	while (pos < TEXT.size()){
		while (pos < TEXT.size() && TEXT[pos] != '(') pos++;
		if (pos < TEXT.size()) detectSupport();
	}

	if (maxv >= 2){
		LOG << "Bootstrap-like support values (0-100) detected.\n";
		maxv = 100;
		defaultv = minv = 0;
	}
	else if (minv >= 0.3) {
		LOG << "Local-probability-like support values (0.333-1) detected.\n";
		maxv = 1;
		defaultv = minv = 0.333;
	}
	else {
		LOG << "Likelihood-like support values (0-1) detected.\n";
		maxv = 1;
		defaultv = minv = 0;
	}

	TEXT = "";
	pos = 0;
}

int main(int argc, char** argv){
	ARG.setProgramName("wastral", "Weighted ASTRAL (all versions included) \n *** You may use --mode to switch from hybrid weighting to other weighting schemes. ***");
	ARG.addStringArg('a', "mapping", "", "A list of gene name to taxon name maps, each line contains one gene name followed by one taxon name separated by a space or tab");
	ARG.addDoubleArg('x', "max", -1, "Max possible support value in weight scale (auto-detect by default)");
	ARG.addDoubleArg('n', "min", -1, "Min possible support value in weight scale (auto-detect by default)");
	ARG.addDoubleArg('d', "default", 0, "Default support value when weight not provided");
	ARG.addFlag('S', "bootstrap", "Bootstrap support value mode (default mode, `-x 100 -n 0 -d 0`)", [&](){
		ARG.getDoubleArg("max") = 100; ARG.getDoubleArg("min") = 0; ARG.getDoubleArg("default") = 0;
	}, true);
	ARG.addFlag('L', "lrt", "Likelihood (alrt) support value mode (`-x 1 -n 0 -d 0`)", [&](){
		ARG.getDoubleArg("max") = 1; ARG.getDoubleArg("min") = 0; ARG.getDoubleArg("default") = 0;
	}, true);
	ARG.addFlag('B', "bayes", "Probability (abayes) support value mode (`-x 1 -n 0.333 -d 0.333`)", [&](){
		ARG.getDoubleArg("max") = 1; ARG.getDoubleArg("min") = 0.333; ARG.getDoubleArg("default") = 0.333;
	}, true);
	ARG.addDoubleArg('w', "weight", -1, "Weight factor of total terminal branch lengths");
	ARG.addIntArg(0, "mode", 1, "1: hybrid weighting, 2: support only, 3: length only, 4: unweighted");
	ARG.addStringArg(0, "treeweights", "", "A file containing a list of gene tree weights, space/tab/new-line separated (default: uniform weights)");

	int dupType = 1;
	string mappingFile;
	meta.initialize(argc, argv);
	mappingFile = ARG.getStringArg("mapping");
	maxv = ARG.getDoubleArg("max");
	minv = ARG.getDoubleArg("min");
	defaultv = ARG.getDoubleArg("default");
	lengthFactor = ARG.getDoubleArg("weight");
	
	if (maxv < 0 || minv < 0) detectSupportType(ARG.getStringArg("input"));

	int mode = ARG.getIntArg("mode");
	useSupport = (mode == 1 || mode == 2);
	useLength = (mode == 1 || mode == 3);

	for (int i = 0; i < meta.nThreads; i++){
		tripInit.nodes.emplace_back();
		tripInit.leafParent.emplace_back();
	}
	if (ARG.getStringArg("root") != ""){
		string s = ARG.getStringArg("root");
		if (name2id.count(s) == 0){
			name2id[s] = names.size();
			names.push_back(s);
			for (int i = 0; i < meta.nThreads; i++){
				tripInit.leafParent[i].emplace_back();
			}
		}
	}
	readInputTrees(ARG.getStringArg("input"), mappingFile, ARG.getStringArg("treeweights"));
	
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
	LOG << "#Genetrees: " << K << endl;
	
	score_t score = meta.run().first;
	//LOG << "Score: " << to_string(score) << endl;
	fprintf(stderr, "Score: %.10lg\n", (double) score);
	return 0;
}

#define DRIVER_VERSION "8"

/* CHANGE LOG
 * 8: bug fix for huge input
 * 7: Adding support for genic branch length
 * 6: Adding compiler flag for coalescent unit
 * 5: Integrating CASTLES-II
 * 4: Adding rooting option
 * 3: Parse gene tree branch lengths
 * 2: Officialize program name
 * 1: Use genetreewithbinaryweight.hpp instead of genetree.hpp 
 */

#include<iostream>
#include<fstream>
#include<unordered_map>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<string>

using namespace std;

#define CASTLES

//#define LARGE_DATA
#ifdef LARGE_DATA
typedef __int128 score_t;
typedef long double length_t;

string to_string(const __int128 x){
	return to_string((double) x);
}

ostream& operator<<(ostream& cout, __int128 x){
	return cout << to_string(x);
}
#else
typedef long long score_t;
typedef double length_t;
#endif

double from_string(const string s){
	return stold(s);
}

#include "argparser.hpp"
#include "genetreewithbinaryweight.hpp"
#include "treeutils.hpp"
#include "castles.hpp"
#include "algorithms.hpp"

MetaAlgorithm meta;
TripartitionInitializer &tripInit = meta.tripInit;

string TEXT;
long long pos = 0;
int K = 0;
int part = 0, iBatch = 0;
vector<string> &names = meta.names;
unordered_map<string, int> &name2id = meta.name2id;
unordered_map<string, string> &leafname_mapping = meta.leafname_mapping;

int MAPPING(int begin, int end){
	string s;
	for (int i = begin; i < end && TEXT[i] != ':'; i++){
		if (TEXT[i] != '\"' && TEXT[i] != '\'') s += TEXT[i];
	}
	if (leafname_mapping.count(s)) s = leafname_mapping[s];
	if (name2id.count(s) == 0){
		name2id[s] = names.size();
		names.push_back(s);
		for (int i = 0; i < meta.nThreads; i++){
			tripInit.leafParent[i].emplace_back();
		}
	}
	return name2id[s];
}

double LENGTH(int begin, int end){
	int i = begin;
	while (i < end && TEXT[i] != ':') i++;
	if (i == end) return 0;
	else return from_string(TEXT.substr(i + 1, end - i - 1));
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
		tripInit.nodes[part][cur].isGhostBranch = true;
		tripInit.nodes[part][cur].length = 0;
		while (TEXT[pos] != ')'){
			int left = lst[rand() % lst.size()];
			int up = tripInit.nodes[part].size();
			tripInit.nodes[part].emplace_back();
			lst.push_back(up);
			if (cur == left) cur = up;
			tripInit.nodes[part][up].isGhostBranch = true;
			tripInit.nodes[part][up].length = 0;
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
		tripInit.nodes[part][cur].isGhostBranch = false;
		tripInit.nodes[part][cur].length = LENGTH(i, pos);
	} 
	else {
		int i = pos;
		while (TEXT[pos] != ')' && TEXT[pos] != ',') pos++;
		tripInit.leafParent[part][MAPPING(i, pos)].push_back(cur);
		tripInit.nodes[part][cur].isGhostBranch = false;
		tripInit.nodes[part][cur].length = LENGTH(i, pos);
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

void examplePrintSubtreeWithSupport(shared_ptr<AnnotatedTree::Node> node){
	if (node->isLeaf()){
		cout << node->taxonName();
	}
	else {
		cout << "(";
		examplePrintSubtreeWithSupport(node->leftChild());
		cout << ",";
		examplePrintSubtreeWithSupport(node->rightChild());
		cout << ")'{support:" << node->support() << ",length:" << node->length()
			<< ",LR_SO:{quartetCnt:" << node->annotation().ab_cd.quartetCnt << ",sumInternal:" <<  node->annotation().ab_cd.sumInternalLength
			<< ",sumL:" <<  node->annotation().ab_cd.sumLengthD << ",sumR:" <<  node->annotation().ab_cd.sumLengthC << ",sumS:" <<  node->annotation().ab_cd.sumLengthB << ",sumO:" <<  node->annotation().ab_cd.sumLengthA << "}"
			<< ",LS_RO:{quartetCnt:" << node->annotation().ac_bd.quartetCnt << ",sumInternal:" <<  node->annotation().ac_bd.sumInternalLength
			<< ",sumL:" <<  node->annotation().ac_bd.sumLengthD << ",sumR:" <<  node->annotation().ac_bd.sumLengthB << ",sumS:" <<  node->annotation().ac_bd.sumLengthC << ",sumO:" <<  node->annotation().ac_bd.sumLengthA << "}"
			<< ",LO_RS:{quartetCnt:" << node->annotation().ad_bc.quartetCnt << ",sumInternal:" <<  node->annotation().ad_bc.sumInternalLength
			<< ",sumL:" <<  node->annotation().ad_bc.sumLengthB << ",sumR:" <<  node->annotation().ad_bc.sumLengthD << ",sumS:" <<  node->annotation().ad_bc.sumLengthC << ",sumO:" <<  node->annotation().ad_bc.sumLengthA << "}}'";
	}
}

string sampleGenetree(int p, int inode, const unordered_set<int> &outgroup){
	string str;
	if (tripInit.nodes[p][inode].small == -1){
		if (outgroup.count(inode)) str = "outgroup";
		else str = "ingroup";
	}
	else str = string("(") + sampleGenetree(p, tripInit.nodes[p][inode].large, outgroup) + "," + sampleGenetree(p, tripInit.nodes[p][inode].small, outgroup) + ")";
	return str + ":" + to_string(tripInit.nodes[p][inode].length);
}

pair<double, int> outgroupLength(const Tree& tr, string label = "outgroup"){
	vector<int> rootChildren;
	int outId = -1;
	for (int i = 0; i < tr.nodes.size(); i++){
		if (tr.nodes[i].label == label) outId = i;
		if (tr.nodes[i].parent == tr.root) rootChildren.push_back(i);
	}
	if (outId == -1 || rootChildren.size() < 2) return {0.0, 0};
	if (rootChildren.size() > 2 || tr.nodes[outId].parent != tr.root) return {tr.nodes[outId].length, 1};
	return {tr.nodes[rootChildren[0]].length + tr.nodes[rootChildren[1]].length, 1};
}

double computeOutgroupLength(){
	double totalLength = 0;
	int cnt = 0;
	for (int p = 0; p < meta.nThreads; p++){
		unordered_set<int> outgroup;
		for (int e: tripInit.leafParent[p][0]) outgroup.insert(e);
		for (int e = 0; e < tripInit.nodes[p].size(); e++){
			if (tripInit.nodes[p][e].up == -1) {
				pair<double, int> treeLength = outgroupLength(sampleGenetree(p, e, outgroup) + ";");
				totalLength += treeLength.first;
				cnt += treeLength.second;
			}
		}
	}
	return totalLength / cnt;
}

int main(int argc, char** argv){
	ARG.setProgramName("astral4", "Accurate Species TRee ALgorithm IV (ASTRAL-IV)\n*** NOW with integrated CASTLES-2 ***");
	ARG.addDoubleArg(0, "genelength", 1000, "Average gene sequence length");
	ARG.addStringArg(0, "length", "SULength", "SULength: substitution-per-site unit; CULength: coalescent unit");
	ARG.addIntArg(0, "sulengthtype", 1, "1: species tree branch length; 2: mean gene tree branch length");
	meta.initialize(argc, argv);
	ARG.getStringArg("annotation") = "localPP";
	
	string mappingFile = ARG.getStringArg("mapping");

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
	readInputTrees(ARG.getStringArg("input"), mappingFile);
	
	LOG << "#Genetrees: " << K << endl;
	
	ARG.getDoubleArg("outgrouplength") = computeOutgroupLength();
	ARG.getIntArg("numgenetrees") = K;

	score_t score = meta.run().first;
	LOG << "Score: " << score << endl;
	return 0;
}

#define DRIVER_VERSION "3"

/* CHANGE LOG
 * 3: Add support for polytomies
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
#include "multitree.hpp"
#include "algorithms.hpp"

const bool VERBOSE = true;

MetaAlgorithm meta;
TripartitionInitializer &tripInit = meta.tripInit;
vector<TripartitionInitializer> &batchInit = meta.batchInit;

unordered_map<string, string> leafname_mapping;
string TEXT;
int duploss = 0;
int nodecnt = 0;
int pos = 0;
int K = 0;
vector<string> &id2name = meta.names;
string rootNtagTrees;
unordered_map<string, int> &name2id = meta.name2id;
bool resolvePolytomies = false, rootNtag = false;
	
class DynamicBitset{
	int size = 0;
	vector<uint64_t> vec;
	
public:
	DynamicBitset(){}
	DynamicBitset(int sz): size(sz), vec((sz + 63) / 64){}
	
	void set(int i){
		if (i >= size){
			size = i + 1;
			if ((size + 63) / 64 > vec.size()){
				vec.resize((size + 63) / 64);
			}
		}
		vec[i / 64] |= (1LL << (i % 64));
	}
	
	DynamicBitset operator|(const DynamicBitset &b) const{
		if (size < b.size) return b | *this;
		DynamicBitset res(size);
		for (int i = 0; i < b.vec.size(); i++){
			res.vec[i] = vec[i] | b.vec[i];
		}
		for (int i = b.vec.size(); i < vec.size(); i++){
			res.vec[i] = vec[i];
		}
		return res;
	}
	
	DynamicBitset operator&(const DynamicBitset &b) const{
		if (size < b.size) return b & *this;
		DynamicBitset res(b.size);
		for (int i = 0; i < b.vec.size(); i++){
			res.vec[i] = vec[i] & b.vec[i];
		}
		return res;
	}
	
	DynamicBitset operator^(const DynamicBitset &b) const{
		if (size < b.size) return b ^ *this;
		DynamicBitset res(size);
		for (int i = 0; i < b.vec.size(); i++){
			res.vec[i] = vec[i] ^ b.vec[i];
		}
		for (int i = b.vec.size(); i < vec.size(); i++){
			res.vec[i] = vec[i];
		}
		return res;
	}
	
	DynamicBitset operator-(const DynamicBitset &b) const{
		DynamicBitset res(size);
		for (int i = 0; i < vec.size(); i++){
			if (i < b.vec.size()) res.vec[i] = vec[i] & ~b.vec[i];
			else res.vec[i] = vec[i];
		}
		return res;
	}
	
	bool operator==(const DynamicBitset &b) const{
		if (size < b.size) return b == *this;
		for (int i = 0; i < b.vec.size(); i++){
			if (vec[i] != b.vec[i]) return false;
		}
		for (int i = b.vec.size(); i < vec.size(); i++){
			if (vec[i] != 0) return false;
		}
		return true;
	}
	
	bool operator!=(const DynamicBitset &b) const{
		return !(*this == b);
	}
	
	bool isDisjointTo(const DynamicBitset &b) const{
		if (size < b.size) return b.isDisjointTo(*this);
		for (int i = 0; i < b.vec.size(); i++){
			if ((vec[i] & b.vec[i]) != 0) return false;
		}
		return true;
	}
	
	vector<int> setBits() const{
		vector<int> res;
		for (int i = 0; i < vec.size(); i++){
			for (int j = 0; j < 64; j++){
				if (vec[i] & (1LL << j)) res.push_back(64 * i + j);
			}
		}
		return res;
	}
};

class GenetreeAnnotator{
private:
	struct Node{
		int leftChildId = -1, rightChildId = -1;
		DynamicBitset label;
		bool isDuplication = false;
		bool isPolytomy = false;
		bool isLeaf = false;
		int score = -1;
		int leafId = -1;
		string name;
	};
	
	vector<Node> node;
	
	tuple<int, int, int> createSubtree(const unordered_map<long long, string> &leafname, const unordered_map<long long, tuple<long long, long long, bool> > &children, const long long cur, vector<int> &rootId){
		if (children.count(cur) == 0){
			int curId = node.size();
			node.emplace_back();
			node[curId].isLeaf = true;
			node[curId].score = 0;
			node[curId].name = leafname.at(cur);
			string name = leafname_mapping.count(leafname.at(cur)) ? leafname_mapping[leafname.at(cur)] : leafname.at(cur);
			if (name2id.count(name) == 0){
				name2id[name] = id2name.size();
				id2name.push_back(name);
				for (int p = 0; p < tripInit.nodes.size(); p++){
					tripInit.leafParent[p].emplace_back();
				}
			}
			node[curId].leafId = name2id[name];
			node[curId].label.set(name2id[name]);
			return make_tuple(curId, -1, -1);
		}
		
		tuple<int, int, int> left = createSubtree(leafname, children, get<0>(children.at(cur)), rootId), right = createSubtree(leafname, children, get<1>(children.at(cur)), rootId);
		int cur0 = node.size();
		node.emplace_back();
		int cur1 = node.size();
		node.emplace_back();
		int cur2 = node.size();
		node.emplace_back();
		node[cur0].leftChildId = get<0>(left);
		node[cur0].rightChildId = get<0>(right);
		
		node[cur1].leftChildId = get<0>(left);
		node[cur2].leftChildId = get<0>(right);
		if (get<1>(left) != -1) node[get<1>(left)].rightChildId = cur2;
		if (get<2>(left) != -1) node[get<2>(left)].rightChildId = cur2;
		if (get<1>(right) != -1) node[get<1>(right)].rightChildId = cur1;
		if (get<2>(right) != -1) node[get<2>(right)].rightChildId = cur1;
		
		node[cur0].isPolytomy = node[cur1].isPolytomy = node[cur2].isPolytomy = get<2>(children.at(cur));

		int root1 = node.size();
		node.emplace_back();
		node[root1].leftChildId = get<0>(left);
		node[root1].rightChildId = cur2;
		rootId.push_back(root1);
		int root2 = node.size();
		node.emplace_back();
		node[root2].leftChildId = get<0>(right);
		node[root2].rightChildId = cur1;
		rootId.push_back(root2);
		
		return make_tuple(cur0, cur1, cur2);
	}
	
	int scoreSubtree(int cur){
		if (node[cur].score != -1) return node[cur].score;
		node[cur].score = scoreSubtree(node[cur].leftChildId) + scoreSubtree(node[cur].rightChildId);
		node[cur].label = node[node[cur].leftChildId].label | node[node[cur].rightChildId].label;
		if (node[cur].isPolytomy) node[cur].isDuplication = true;
		else if (!node[node[cur].leftChildId].label.isDisjointTo(node[node[cur].rightChildId].label)){
			node[cur].score++;
			node[cur].isDuplication = true;
			if (node[cur].label != node[node[cur].leftChildId].label) node[cur].score++;
			if (node[cur].label != node[node[cur].rightChildId].label) node[cur].score++;
		}
		return node[cur].score;
	}
	
public:
	const vector<string> &leafnames() const{
		return id2name;
	}
	
	int annotateTree(const unordered_map<long long, string> &leafname, const unordered_map<long long, tuple<long long, long long, bool> > &children, const long long root){
		vector<int> rootId;
		tuple<int, int, int> left = createSubtree(leafname, children, get<0>(children.at(root)), rootId), right = createSubtree(leafname, children, get<1>(children.at(root)), rootId);
		if (get<1>(left) != -1) node[get<1>(left)].rightChildId = get<0>(right);
		if (get<2>(left) != -1) node[get<2>(left)].rightChildId = get<0>(right);
		if (get<1>(right) != -1) node[get<1>(right)].rightChildId = get<0>(left);
		if (get<2>(right) != -1) node[get<2>(right)].rightChildId = get<0>(left);
		
		int curId = node.size();
		node.emplace_back();
		node[curId].leftChildId = get<0>(left);
		node[curId].rightChildId = get<0>(right);
		rootId.push_back(curId);
		
		int bestscore = 999999, bestroot = -1, bestcnt = 0;
		for (int root: rootId){
			int score = scoreSubtree(root);
			if (score == bestscore){
				bestcnt++;
				if (rand() % bestcnt == 0) bestroot = root;
			}
			if (score < bestscore){
				bestscore = score;
				bestroot = root;
				bestcnt = 1;
			}
		}
		duploss += bestscore;
		return bestroot;
	}
	
	int buildTree(int cur, int p) const{
		int w = tripInit.nodes[p].size();
		tripInit.nodes[p].emplace_back();
		
		if (node[cur].isLeaf){
			tripInit.leafParent[p][node[cur].leafId].push_back(w);
			//cerr << id2name[node[cur].leafId];
			return w;
		}
		int left = node[cur].leftChildId, right = node[cur].rightChildId;
		int large, small;
		if (node[left].label.setBits().size() > node[right].label.setBits().size()){
			large = left;
			small = right;
		}
		else {
			large = right;
			small = left;
		}
		//cerr << "(";
		int u = buildTree(small, p);
		//cerr << ",";
		int v = buildTree(large, p);
		//cerr << ")";
		tripInit.nodes[p][w].small = u;
		tripInit.nodes[p][w].large = v;
		tripInit.nodes[p][u].up = w;
		tripInit.nodes[p][v].up = w;
		if (node[cur].isDuplication){
			tripInit.nodes[p][w].dup = true;
			for (int i: (node[cur].label - node[large].label).setBits()){
				tripInit.leafParent[p][i].push_back(w);
			} //cerr << "+";
		}
		else {
			tripInit.nodes[p][w].dup = false;
		}
		return w;
	}

	string printTree(int cur, bool isRoot = true) const{
		if (node[cur].isLeaf){
			if (isRoot) return node[cur].name + ";";
			return node[cur].name;
		}
		int left = node[cur].leftChildId, right = node[cur].rightChildId;
		string res = string("(") + printTree(left, false) + "," + printTree(right, false) + ")";
		if (node[cur].isDuplication) res += "D";
		if (isRoot) res += ";";
		return res;
	}
};

string MAPPING(int begin, int end){
	string s;
	for (int i = begin; i < end && TEXT[i] != ':'; i++){
		if (TEXT[i] != '\"' && TEXT[i] != '\'') s += TEXT[i];
	}
	if (leafname_mapping.count(s)) return leafname_mapping[s];
	else return s;
}

string GET_NAME(int begin, int end){
	string s;
	for (int i = begin; i < end && TEXT[i] != ':'; i++){
		if (TEXT[i] != '\"' && TEXT[i] != '\'') s += TEXT[i];
	}
	return s;
}

long long parse(unordered_map<long long, string> &leafname, unordered_map<long long, tuple<long long, long long, bool> > &children, bool isRoot = false){
	int i = pos;
	long long cur;
	while (TEXT[pos] != '(' && TEXT[pos] != ',' && TEXT[pos] != ')') pos++;
	if (TEXT[pos] != '(') {
		cur = nodecnt++;
		leafname[cur] = GET_NAME(i, pos);
		return cur;
	}
	else {
		pos++; // (
		cur = parse(leafname, children);
		while (TEXT[pos] != ',') pos++;
		vector<long long> lst, plst;
		lst.push_back(cur);
		while (TEXT[pos] != ')'){
			pos++; // ,
			long long temp = lst[rand() % lst.size()];
			long long left = nodecnt++, right = parse(leafname, children);
			if (leafname.count(temp)) {
				leafname[left] = leafname[temp];
				leafname.erase(temp);
			}
			else children[left] = children[temp];
			children[temp] = make_tuple<long long, long long, bool>((long long) left, (long long) right, true);
			plst.push_back(temp);
			lst.push_back(left);
			lst.push_back(right);
		}
		if ( (lst.size() == 5 && isRoot) || lst.size() == 3){
			for (long long temp: plst) get<2>(children[temp]) = false;
		}
		while (TEXT[pos] != ')') pos++;
		pos++; // )
		while (TEXT[pos] != ',' && TEXT[pos] != ')' && TEXT[pos] != ';') pos++;
		return cur;
	}
}

string convert2string(const unordered_map<long long, string> &leafname, const unordered_map<long long, tuple<long long, long long, bool> > &children, long long node){
	if (leafname.count(node)) return leafname.at(node);
	const tuple<long long, long long, bool> &e = children.at(node);
	return string("(") + convert2string(leafname, children, get<0>(e)) + "," + convert2string(leafname, children, get<1>(e)) + ")"
		+ (get<2>(children.at(node)) ? "P" : "");
}

void annotate(string input, string mapping){
	bool hasPolytomy = false;
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
			int leafCnt = 1, internalCnt = 0;
			for (int i = pos; i < TEXT.size() && TEXT[i] != ';'; i++){
				if (TEXT[i] == '(') internalCnt++;
				if (TEXT[i] == ',') leafCnt++;
			}
			if (internalCnt < leafCnt - 2 && !resolvePolytomies && !hasPolytomy) {
				cerr << "Warning: Non-binary input tree(s) detected! ASTRAL-Pro will treat all polytomies as duplication events. Please ignore this warning if all polytomies are parents of leaves (eg. output of Fasttree).\n";
				hasPolytomy = true;
			}
			
			unordered_map<long long, string> leafname;
			unordered_map<long long, tuple<long long, long long, bool> > children;
			long long root = parse(leafname, children, true);
			GenetreeAnnotator ga;
			//cerr << convert2string(leafname, children, root) << endl;
			int iroot = ga.annotateTree(leafname, children, root);
			if (rootNtag) rootNtagTrees += ga.printTree(iroot) + "\n";
			else ga.buildTree(iroot, K % tripInit.nodes.size());
			K++;
			if (VERBOSE && (K & 511) == 0) cerr << "Read " << K << " genetrees and found " << id2name.size() << " species.\n";
		}
	}
}

string HELP_TEXT = R"V0G0N(-a  a list of gene name to taxon name maps, each line contains one gene name followed by one taxon name separated by a space or tab 
-e  0(default): exit when input contains polytomies; 1: resolve polytomies (no guarentee on accuracy)
inputGeneTrees: the path to a file containing all gene trees in Newick format
)V0G0N";

int main(int argc, char** argv){
	ARG.setProgramName("astral-pro", "ASTRAL for PaRalogs and Orthologs");
	ARG.addStringArg('a', "mapping", "", "A list of gene name to taxon name maps, each line contains one gene name followed by one taxon name separated by a space or tab", true);
	ARG.addIntArg('e', "exit", 0, "0: print warning when input contains polytomies; 1: resolving polytomies; 2: printing rooted and tagged gene trees and exit");
	ARG.addFlag('E', "noexit", "No warning when input contains polytomies (`-e 1`)", [&](){
			ARG.getIntArg("exit") = 1;
	});
	ARG.addFlag('T', "tagging", "Just printing rooted and tagged gene trees (`-e 2`)", [&](){
			ARG.getIntArg("exit") = 2;
	});
	string mappingFile;
	meta.initialize(argc, argv, " -a taxonNameMaps", HELP_TEXT);
	mappingFile = ARG.getStringArg("mapping");
	resolvePolytomies = (ARG.getIntArg("exit") != 0);
	rootNtag = (ARG.getIntArg("exit") == 2);
	/*
	for (int i = 1; i < argc; i += 2){
		if (strcmp(argv[i], "-y") == 0) {i--; continue;}
		if (strcmp(argv[i], "-e") == 0) resolvePolytomies = stoi(argv[i + 1]);
		if (strcmp(argv[i], "-a") == 0) mappingFile = argv[i + 1];
	}
	*/

	for (int i = 0; i < meta.nThread2; i++){
		tripInit.nodes.emplace_back();
		tripInit.leafParent.emplace_back();
	}
	annotate(ARG.getStringArg("input"), mappingFile);
	
	cerr << "#Genetrees: " << K << endl;
	cerr << "#Duploss: " << duploss << endl;
	
	if (rootNtag){
		if (ARG.getStringArg("output") == "<standard output>") cout << rootNtagTrees;
		else ofstream(ARG.getStringArg("output")) << rootNtagTrees;
		exit(0);
	}

	score_t score = meta.run().first;
	cerr << "#EqQuartets: " << NUM_EQ_CLASSES << endl;
	cerr << "Score: " << score << endl;
	return 0;
}

#define DRIVER_VERSION "6"

/* CHANGE LOG
 * 6: Adding compiler flag for coalescent unit
 * 5: Integrate CASTLES-Pro
 * 4: Add support for CASTLES
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
#include "multitree.hpp"
#include "treeutils.hpp"
#include "castles.hpp"
#include "algorithms.hpp"

const bool VERBOSE = true;

MetaAlgorithm meta;
TripartitionInitializer &tripInit = meta.tripInit;

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
		length_t length = 0;
	};
	
	vector<Node> node;
	
	tuple<int, int, int> createSubtree(const unordered_map<long long, string> &leafname, const unordered_map<long long, tuple<long long, long long, bool> > &children, 
			const unordered_map<long long, length_t> &edgelength, const long long cur, vector<int> &rootId){
		if (children.count(cur) == 0){
			int curId = node.size();
			node.emplace_back();
			node[curId].length = edgelength.at(cur);
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
		
		tuple<int, int, int> left = createSubtree(leafname, children, edgelength, get<0>(children.at(cur)), rootId), right = createSubtree(leafname, children, edgelength, get<1>(children.at(cur)), rootId);
		int cur0 = node.size();
		node.emplace_back();
		int cur1 = node.size();
		node.emplace_back();
		int cur2 = node.size();
		node.emplace_back();
		
		node[cur0].length = edgelength.at(cur);
		node[cur1].length = edgelength.at(get<1>(children.at(cur)));
		node[cur2].length = edgelength.at(get<0>(children.at(cur)));

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
	
	int annotateTree(const unordered_map<long long, string> &leafname, const unordered_map<long long, tuple<long long, long long, bool> > &children, const unordered_map<long long, length_t> &edgelength, const long long root){
		vector<int> rootId;
		tuple<int, int, int> left = createSubtree(leafname, children, edgelength, get<0>(children.at(root)), rootId), right = createSubtree(leafname, children, edgelength, get<1>(children.at(root)), rootId);
	
		length_t length = node[get<0>(left)].length + node[get<0>(right)].length;
		node[get<0>(left)].length = length;
		node[get<0>(right)].length = length;

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
	
	int buildTree(int cur, int p, int depth = 0) const{
		int w = tripInit.nodes[p].size();
		tripInit.nodes[p].emplace_back();
		
		if (node[cur].isLeaf) tripInit.leafParent[p][node[cur].leafId].push_back(w);
		else{
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
			int u = buildTree(small, p, depth + 1);
			int v = buildTree(large, p, depth + 1);
			tripInit.nodes[p][w].small = u;
			tripInit.nodes[p][w].large = v;
			tripInit.nodes[p][u].up = w;
			tripInit.nodes[p][v].up = w;
			if (node[cur].isDuplication){
				tripInit.nodes[p][w].dup = true;
				for (int i: (node[cur].label - node[large].label).setBits()){
					tripInit.leafParent[p][i].push_back(w);
				}
			}
			else {
				tripInit.nodes[p][w].dup = false;
			}
		}
		tripInit.nodes[p][w].length = ((depth == 1) ? node[cur].length / 2 : node[cur].length);
		return w;
	}

	string printTree(int cur, int depth = 0) const{
		string res;
		if (node[cur].isLeaf) res = node[cur].name;
		else {
			int left = node[cur].leftChildId, right = node[cur].rightChildId;
			res = string("(") + printTree(left, depth + 1) + "," + printTree(right, depth + 1) + ")";
			if (node[cur].isDuplication) res += "D";
		}
		if (depth == 0) res += ";";
		else if (depth == 1) res += string(":") + to_string(node[cur].length / 2);
		else res += string(":") + to_string(node[cur].length);
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

double LENGTH(int begin, int end){
	int i = begin;
	while (i < end && TEXT[i] != ':') i++;
	if (i == end) return 0;
	else return from_string(TEXT.substr(i + 1, end - i - 1));
}

long long parse(unordered_map<long long, string> &leafname, unordered_map<long long, tuple<long long, long long, bool> > &children, unordered_map<long long, length_t> &edgelength, bool isRoot = false){
	int i = pos;
	long long cur;
	while (TEXT[pos] != '(' && TEXT[pos] != ',' && TEXT[pos] != ')') pos++;
	if (TEXT[pos] != '(') {
		cur = nodecnt++;
		leafname[cur] = GET_NAME(i, pos);
		edgelength[cur] = LENGTH(i, pos);
		return cur;
	}
	else {
		pos++;
		cur = parse(leafname, children, edgelength);
		while (TEXT[pos] != ',') pos++;
		vector<long long> lst, plst;
		lst.push_back(cur);
		while (TEXT[pos] != ')'){
			pos++;
			long long temp = lst[rand() % lst.size()];
			long long left = nodecnt++, right = parse(leafname, children, edgelength);
			if (leafname.count(temp)) {
				leafname[left] = leafname[temp];
				leafname.erase(temp);
			}
			else children[left] = children[temp];
			edgelength[left] = edgelength[temp];
			children[temp] = make_tuple<long long, long long, bool>((long long) left, (long long) right, true);
			plst.push_back(temp);
			lst.push_back(left);
			lst.push_back(right);
		}
		if ( (lst.size() == 5 && isRoot) || lst.size() == 3){
			for (long long temp: plst) get<2>(children[temp]) = false;
		}
		while (TEXT[pos] != ')') pos++;
		pos++;
		i = pos;
		while (TEXT[pos] != ',' && TEXT[pos] != ')' && TEXT[pos] != ';') pos++;
		edgelength[cur] = LENGTH(i, pos);
		return cur;
	}
}

string convert2string(const unordered_map<long long, string> &leafname, const unordered_map<long long, tuple<long long, long long, bool> > &children, const unordered_map<long long, length_t> &edgelength, long long node){
	if (leafname.count(node)) return leafname.at(node) + ":" + to_string(edgelength.at(node));
	const tuple<long long, long long, bool> &e = children.at(node);
	return string("(") + convert2string(leafname, children, edgelength, get<0>(e)) + "," + convert2string(leafname, children, edgelength, get<1>(e)) + ")"
		+ (get<2>(children.at(node)) ? "P" : "") + ":" + to_string(edgelength.at(node));
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
				LOG << "Warning: Non-binary input tree(s) detected! ASTRAL-Pro will treat all polytomies as duplication events. Please ignore this warning if all polytomies are parents of leaves (eg. output of Fasttree).\n";
				hasPolytomy = true;
			}
			
			unordered_map<long long, string> leafname;
			unordered_map<long long, tuple<long long, long long, bool> > children;
			unordered_map<long long, length_t> edgelength;
			long long root = parse(leafname, children, edgelength, true);
			GenetreeAnnotator ga;
			// LOG << convert2string(leafname, children, edgelength, root) << endl;
			int iroot = ga.annotateTree(leafname, children, edgelength, root);
			// LOG << ga.printTree(iroot) << "\n";
			if (rootNtag) rootNtagTrees += ga.printTree(iroot) + "\n";
			else ga.buildTree(iroot, K % tripInit.nodes.size());
			K++;
			if (VERBOSE && (K & 511) == 0) LOG << "Read " << K << " genetrees and found " << id2name.size() << " species.\n";
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
	else if (tripInit.nodes[p][inode].dup){
		if (rand() & 1) return sampleGenetree(p, tripInit.nodes[p][inode].large, outgroup);
		else return sampleGenetree(p, tripInit.nodes[p][inode].small, outgroup);
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
	ARG.setProgramName("astral-pro3", "ASTRAL for PaRalogs and Orthologs III (ASTRAL-Pro3)\n*** NOW with integrated CASTLES-Pro ***");
	ARG.addDoubleArg(0, "genelength", 1000, "Average gene sequence length");
	ARG.addStringArg(0, "length", "SULength", "SULength: substitution-per-site unit; CULength: coalescent unit");
	ARG.addIntArg('e', "exit", 0, "0: print warning when input contains polytomies; 1: resolving polytomies; 2: printing rooted and tagged gene trees and exit");
	ARG.addFlag('E', "noexit", "No warning when input contains polytomies (`-e 1`)", [&](){
			ARG.getIntArg("exit") = 1;
	});
	ARG.addFlag('T', "tagging", "Just printing rooted and tagged gene trees (`-e 2`)", [&](){
			ARG.getIntArg("exit") = 2;
	});
	meta.initialize(argc, argv);
	ARG.getStringArg("annotation") = "localPP";

	string mappingFile = ARG.getStringArg("mapping");
	resolvePolytomies = (ARG.getIntArg("exit") != 0);
	rootNtag = (ARG.getIntArg("exit") == 2);

	for (int i = 0; i < meta.nThreads; i++){
		tripInit.nodes.emplace_back();
		tripInit.leafParent.emplace_back();
	}
	if (ARG.getStringArg("root") != ""){
		string s = ARG.getStringArg("root");
		if (name2id.count(s) == 0){
			name2id[s] = id2name.size();
			id2name.push_back(s);
			for (int i = 0; i < meta.nThreads; i++){
				tripInit.leafParent[i].emplace_back();
			}
		}
	}
	annotate(ARG.getStringArg("input"), mappingFile);
	
	LOG << "#Genetrees: " << K << endl;
	LOG << "#Duploss: " << duploss << endl;
	
	if (rootNtag){
		if (ARG.getStringArg("output") == "<standard output>") cout << rootNtagTrees;
		else ofstream(ARG.getStringArg("output")) << rootNtagTrees;
		exit(0);
	}

	#ifdef CASTLES
	ARG.getDoubleArg("outgrouplength") = computeOutgroupLength();
	ARG.getIntArg("numgenetrees") = K;
	#endif

	score_t score = meta.run().first;
	LOG << "#EqQuartets: " << NUM_EQ_CLASSES << endl;
	LOG << "Score: " << score << endl;

	/*
	#ifdef CASTLES
	examplePrintSubtreeWithSupport(meta.annotTree->root());
	cout << ";\n";
	#endif
	*/
	return 0;
}

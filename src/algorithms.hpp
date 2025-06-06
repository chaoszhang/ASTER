#define ALG_VERSION "v1.22"

/* CHANGE LOG
 * 1.22: adding back freqQuad.csv
 * 1.21: de-warning and new species tree representation
 * 1.20: updating NNI algorithm
 * 1.19: error message for guide trees and constraint trees
 * 1.18: improving numerical stability for incomplete beta
 * 1.17: adding support for castles
 * 1.16: adding support for annotations for branch lengths
 * 1.15: adding rooting option
 * 1.14: adding more information in -u 2 option
 * 1.13: adding species tree annotation functionality
 * 1.12: significantly refactoring the code for reduced memory
 * 1.11: fixing incomplete class bug in threadpool for some compilers
 * 1.10: fixing constraint/guide trees for rooting on nodes
 * 1.9: deleting obsolete code
 * 1.8: adding seed to pseudorandomness
 * 1.7: adding Grubbs's test based support
 * 1.6: adding -w option for gene tree replications
 * 1.5: improving parallelization
 */

typedef unsigned __int128 hash_t;

#include<iostream>
#include<vector>
#include<array>
#include<utility> 
#include<random>
#include<unordered_set>
#include<unordered_map>
#include<algorithm>
#include<tuple>
#include<string>
#include<thread>
#include<cmath>
#include<memory>
#include<iomanip>

#include "incbeta.c"
#include "argparser.hpp"
#include "threadpool.hpp"
#include "speciestree.hpp"

#ifndef ERROR_TOLERANCE
#define ERROR_TOLERANCE 0
#endif

#if defined(SUPPORT) || defined(G_SUPPORT) || defined(LOCAL_BOOTSTRAP)
#define SOME_FORM_OF_SUPPORT
#endif

using namespace std;

struct LogInfo {
	ostream& fout;
	bool enabled;

	LogInfo(bool enabled = true, ostream& fout = cerr): fout(fout), enabled(enabled) {}

	template<class T> LogInfo& operator<< (const T& v) {
		if (enabled) fout << v;
		return *this;
	}
	
	LogInfo& operator<< (basic_ostream<char, std::char_traits<char> >& (*endl)(basic_ostream<char, std::char_traits<char> >&)) {
		if (enabled) fout << endl;
		return *this;
	}
} LOG;

mt19937_64 RND_GENERATOR;

struct Hasher{
	size_t operator()(const hash_t h) const{
		return (size_t) h;
	}
};

string formatBootstrap(double s){
	int a = floor(s);
	int b = round(10 * (s - a));
	return to_string(a) + "." + to_string(b);
}

struct AnnotatedTree{
	struct Node{
		shared_ptr<Node> lc, rc;
		weak_ptr<Node> p;
		int taxon = -1;
		string name;
		double s = 0, len = 0;

		#ifdef CUSTOMIZED_ANNOTATION
		CustomizedAnnotation annot;
		
		const CustomizedAnnotation& annotation() const{
			return annot;
		}
		#endif

		Node(){}

		bool isLeaf() const{
			return taxon != -1;
		}

		shared_ptr<Node> parent() const{
			return p.lock();
		}

		shared_ptr<Node> leftChild() const{
			return lc;
		}

		shared_ptr<Node> rightChild() const{
			return rc;
		}

		int taxonID() const{
			return taxon;
		}

		string taxonName() const{
			return name;
		}

		double support() const{
			return s;
		}

		double length() const{
			return len;
		}
	};

	shared_ptr<Node> r;

	shared_ptr<Node> root(){
		return r;
	}

	shared_ptr<Node> addRoot(){
		r.reset(new Node());
		return r;
	}

	shared_ptr<Node> addLeft(shared_ptr<Node> node){
		shared_ptr<Node> child(new Node());
		child->p = node;
		node->lc = child;
		return child;
	}

	shared_ptr<Node> addRight(shared_ptr<Node> node){
		shared_ptr<Node> child(new Node());
		child->p = node;
		node->rc = child;
		return child;
	}
};

struct PlacementAlgorithm{
	
	struct Node{
		Node(int leafId = -1): leafId(leafId){}
		
		int parent = -1, heavy = -1, light = -1; //parent 0, heavy 1, light 2
		int leafId = -1, leafCnt = 1;
		score_t sYes = 0, sNo = 0;
		bool placeHeavy = false, placeLight = false; 
	};
	
	vector<hash_t> taxonHash; // sum of hash = 0
	vector<tuple<hash_t, hash_t, score_t> > tripHash;
	vector<int> order;
	vector<Node> nodes;
	int rootNodeId = -1, rootLeafId = -1;
	size_t orderId = 0;
	int rNN;
	
	Tripartition trip;
	const int ROUND_NN = -1;
	ThreadPool &TP;
	
	PlacementAlgorithm(const vector<hash_t> &taxonHash, TripartitionInitializer &tripInit, ThreadPool& TP, int ROUND_NN):
		taxonHash(taxonHash), rNN(ROUND_NN), trip(tripInit), ROUND_NN(ROUND_NN), TP(TP) {}
	
	int& heavy(int v){
		return nodes[v].heavy;
	}
	
	int& light(int v){
		return nodes[v].light;
	}
	
	int& parent(int v){
		return nodes[v].parent;
	}
	
	int& leafId(int v){
		return nodes[v].leafId;
	}
	
	int& leafCnt(int v){
		return nodes[v].leafCnt;
	}
	
	score_t& sYes(int v){
		return nodes[v].sYes;
	}
	
	score_t& sNo(int v){
		return nodes[v].sNo;
	}
	
	bool& placeHeavy(int v){
		return nodes[v].placeHeavy;
	}
	
	bool& placeLight(int v){
		return nodes[v].placeLight;
	}
	
	int makeNode(int v, int u, bool modifyRoot = false){
		int w = nodes.size();
		nodes.emplace_back(-1);
		parent(w) = parent(v);
		if (parent(v) == -1) {
			if (modifyRoot) rootNodeId = w;
		}
		else if (heavy(parent(w)) == v) heavy(parent(w)) = w;
		else light(parent(w)) = w;
		parent(v) = w; parent(u) = w;
		heavy(w) = v; light(w) = u;
		for (; w != -1; w = parent(w)){
			if (leafCnt(heavy(w)) < leafCnt(light(w))){
				int t = heavy(w);
				heavy(w) = light(w);
				light(w) = t;
			}
			leafCnt(w) = leafCnt(heavy(w)) + leafCnt(light(w));
		}
		return nodes.size() - 1;
	}
	
	void addSibling(int v, int i){
		int u = nodes.size();
		nodes.emplace_back(i);
		makeNode(v, u, true);
	}
	
	void forceLeftHeavy(int v){
		if (leafId(v) != -1) leafCnt(v) = 1;
		else {
			if (leafCnt(heavy(v)) < leafCnt(light(v))){
				int t = heavy(v);
				heavy(v) = light(v);
				light(v) = t;
			}
			leafCnt(v) = leafCnt(heavy(v)) + leafCnt(light(v));
		}
	}
	
	void defaultInitializer(){
		if (rootLeafId == -1){
			rootLeafId = order[orderId];
			orderId++;
		}
		if (rootNodeId == -1){
			rootNodeId = nodes.size();
			nodes.emplace_back(order[orderId]);
			orderId++;
		}
		forceLeftHeavy(rootNodeId);
	}
	
	void place(int i){
		addSibling(locateBranchAdding(i), i);
	}
	
	void tripUpdateSet(int tgt, int i){
		Tripartition &t = trip;
		TP.push([=, &t](int part)->score_t{t.updatePart(part, tgt, i); return 0;});
	}

	void tripUpdateGet(int tgt, int i){
		TP.pop();
	}

	void tripScoreSet(){
		Tripartition &t = trip;
		TP.push([=, &t](int part)->score_t{return t.scorePart(part);});
	}

	score_t tripScoreGet(){
		return TP.pop();
	}

	score_t tripScore(){
		tripScoreSet();
		return tripScoreGet();
	}

	void switchSubtreeSet(int v, int src, int tgt){
		if (leafId(v) != -1){
			tripUpdateSet(tgt, leafId(v));
		}
		else{
			switchSubtreeSet(light(v), src, tgt);
			switchSubtreeSet(heavy(v), src, tgt);
		}
	}
	
	void switchSubtreeGet(int v, int src, int tgt){
		if (leafId(v) != -1){
			tripUpdateGet(tgt, leafId(v));
		}
		else{
			switchSubtreeGet(light(v), src, tgt);
			switchSubtreeGet(heavy(v), src, tgt);
		}
	}

	void switchSubtree(int v, int src, int tgt){
		switchSubtreeSet(v, src, tgt);
		switchSubtreeGet(v, src, tgt);
	}

	void scoringPlacementDPSet(int v, int i){
		if (leafId(v) != -1){
			//Wv|-|i
			//W|v|i
			tripUpdateSet(1, leafId(v));
			tripScoreSet();
		}
		else{
			//AC|B|i
			scoringPlacementDPSet(light(v), i);
			//ABC|-|i
			switchSubtreeSet(light(v), 1, 0);
			//BC|A|i
			scoringPlacementDPSet(heavy(v), i);
			//C|A|Bi
			switchSubtreeSet(light(v), 0, 2);
			tripScoreSet();
			//C|Ai|B
			tripUpdateSet(1, i);
			tripScoreSet();
			//Ci|A|B
			tripUpdateSet(0, i);
			tripScoreSet();
			//C|A|Bi
			tripUpdateSet(2, i);
			//C|AB|i
			switchSubtreeSet(light(v), 2, 1);
			tripScoreSet();
		}
	}
		
	void scoringPlacementDPGet(int v, int i){
		if (leafId(v) != -1){
			//Wv|-|i
			sNo(v) = 0;
			//W|v|i
			tripUpdateGet(1, leafId(v));
			sYes(v) = tripScoreGet();
			placeHeavy(v) = false;
			placeLight(v) = false;
		}
		else{
			//AC|B|i
			scoringPlacementDPGet(light(v), i);
			//ABC|-|i
			switchSubtreeGet(light(v), 1, 0);
			//BC|A|i
			scoringPlacementDPGet(heavy(v), i);
			//C|A|Bi
			switchSubtreeGet(light(v), 0, 2);
			score_t up, left, right, sibling;
			right = tripScoreGet();
			//C|Ai|B
			tripUpdateGet(1, i);
			left = tripScoreGet();
			//Ci|A|B
			tripUpdateGet(0, i);
			up = tripScoreGet();
			//C|A|Bi
			tripUpdateGet(2, i);
			//C|AB|i
			switchSubtreeGet(light(v), 2, 1);
			sibling = tripScoreGet();
			sNo(v) = sNo(light(v)) + sNo(heavy(v)) + up;
			sYes(v) = sNo(v) + sibling;
			placeHeavy(v) = false;
			placeLight(v) = false;
			if (sYes(v) < sYes(light(v)) + sNo(heavy(v)) + right){
				sYes(v) = sYes(light(v)) + sNo(heavy(v)) + right;
				placeHeavy(v) = false;
				placeLight(v) = true;
			}
			if (sYes(v) < sNo(light(v)) + sYes(heavy(v)) + left){
				sYes(v) = sNo(light(v)) + sYes(heavy(v)) + left;
				placeHeavy(v) = true;
				placeLight(v) = false;
			}
		}
	}

	void locateBranchAddingSet(int i){
		tripUpdateSet(2, i);
		tripUpdateSet(0, rootLeafId);
		switchSubtreeSet(rootNodeId, -1, 0);
		scoringPlacementDPSet(rootNodeId, i);
	}
	
	int locateBranchAddingGet(int i){
		//trip.addTotal(i);
		tripUpdateGet(2, i);
		tripUpdateGet(0, rootLeafId);
		switchSubtreeGet(rootNodeId, -1, 0);
		scoringPlacementDPGet(rootNodeId, i);
		int v = rootNodeId;
		while (placeHeavy(v) || placeLight(v)){
			if (placeHeavy(v)) v = heavy(v);
			else v = light(v);
		}
		return v;
	}
	
	int locateBranchAdding(int i){
		locateBranchAddingSet(i);
		return locateBranchAddingGet(i);
	}
	
	int locateBranch(int i){
		locateBranchAddingSet(i);
		tripUpdateSet(-1, i);
		int v = locateBranchAddingGet(i);
		tripUpdateGet(-1, i);
		return v;
	}

	inline score_t nnPerform(int v, int best_case, score_t top_score){
		if (best_case == 1){
			int u = light(v), c = heavy(u), d = light(u), ab = heavy(v);
			parent(ab) = u; heavy(u) = ab;
			parent(c) = u; light(u) = c;
			parent(u) = v; heavy(v) = u;
			parent(d) = v; light(v) = d;
			leafCnt(u) = leafCnt(heavy(u)) + leafCnt(light(u));
			return nnLocal(v, true, top_score);
		}
		if (best_case == 2){
			int u = light(v), c = heavy(u), d = light(u), ab = heavy(v);
			parent(ab) = u; heavy(u) = ab;
			parent(d) = u; light(u) = d;
			parent(u) = v; heavy(v) = u;
			parent(c) = v; light(v) = c;
			leafCnt(u) = leafCnt(heavy(u)) + leafCnt(light(u));
			return nnLocal(v, true, top_score);
		}
		if (best_case == 3){
			int u = heavy(v), a = heavy(u), b = light(u), cd = light(v);
			parent(b) = u; heavy(u) = b;
			parent(cd) = u; light(u) = cd;
			parent(a) = v; heavy(v) = a;
			parent(u) = v; light(v) = u;
			leafCnt(u) = leafCnt(heavy(u)) + leafCnt(light(u));
			return nnLocal(v, true, top_score);
		}
		if (best_case == 4){
			int u = heavy(v), a = heavy(u), b = light(u), cd = light(v);
			parent(a) = u; heavy(u) = a;
			parent(cd) = u; light(u) = cd;
			parent(u) = v; heavy(v) = u;
			parent(b) = v; light(v) = b;
			leafCnt(u) = leafCnt(heavy(u)) + leafCnt(light(u));
			return nnLocal(v, true, top_score);
		}
		return top_score;
	}
	
	score_t nnLocal(int v, bool recurse, score_t e_ab_cd){
		//E|ABCD|-
		if (recurse) {
			if (rNN <= 0) return e_ab_cd;
			rNN--;
		}
		
		switchSubtree(light(v), 1, 2);
		//E|AB|CD
		score_t abe_c_d = 0, cde_a_b = 0;
		if (leafId(light(v)) == -1){
			int c = heavy(light(v)), d = light(light(v));
			switchSubtreeSet(heavy(v), 1, 0);
			switchSubtreeSet(c, 2, 1);
			//ABE|C|D
			tripScoreSet();
			switchSubtreeSet(d, 2, 1);
			//ABE|CD|-

			switchSubtreeGet(heavy(v), 1, 0);
			switchSubtreeGet(c, 2, 1);
			//ABE|C|D
			abe_c_d = tripScoreGet();
			switchSubtreeGet(d, 2, 1);
			//ABE|CD|-
			if (recurse) abe_c_d = nnLocal(light(v), false, abe_c_d);
			switchSubtreeSet(heavy(v), 0, 1);
			switchSubtreeSet(light(v), 1, 2);
			//E|AB|CD

			switchSubtreeGet(heavy(v), 0, 1);
			switchSubtreeGet(light(v), 1, 2);
			//E|AB|CD
		}
		if (leafId(heavy(v)) == -1){
			int b = light(heavy(v));
			//CDE|A|B
			switchSubtreeSet(light(v), 2, 0);
			switchSubtreeSet(b, 1, 2);
			tripScoreSet();
			//CDE|AB|-
			switchSubtreeSet(b, 2, 1);
			
			//CDE|A|B
			switchSubtreeGet(light(v), 2, 0);
			switchSubtreeGet(b, 1, 2);
			cde_a_b = tripScoreGet();
			//CDE|AB|-
			switchSubtreeGet(b, 2, 1);

			if (recurse) cde_a_b = nnLocal(heavy(v), false, cde_a_b);
			//E|AB|CD
			switchSubtree(light(v), 0, 2);
		}
		int best_case = 0;
		score_t best_score = abe_c_d + cde_a_b + e_ab_cd, top_score = e_ab_cd;
		if (leafId(light(v)) == -1){
			int c = heavy(light(v)), d = light(light(v));
			//DE|AB|C
			switchSubtreeSet(d, 2, 0);
			tripScoreSet();
			//E|ABD|C
			switchSubtreeSet(d, 0, 1);
			tripScoreSet();
			//E|ABC|D
			switchSubtreeSet(d, 1, 2);
			switchSubtreeSet(c, 2, 1);
			tripScoreSet();
			//CE|AB|D
			switchSubtreeSet(c, 1, 0);
			tripScoreSet();
			//E|AB|CD
			switchSubtreeSet(c, 0, 2);

			//DE|AB|C
			switchSubtreeGet(d, 2, 0);
			score_t de_ab_c = tripScoreGet();
			//E|ABD|C
			switchSubtreeGet(d, 0, 1);
			score_t e_abd_c = tripScoreGet();
			//E|ABC|D
			switchSubtreeGet(d, 1, 2);
			switchSubtreeGet(c, 2, 1);
			score_t e_abc_d = tripScoreGet();
			//CE|AB|D
			switchSubtreeGet(c, 1, 0);
			score_t ce_ab_d = tripScoreGet();
			if (best_score + ERROR_TOLERANCE < cde_a_b + de_ab_c + e_abc_d) {
				best_score = cde_a_b + de_ab_c + e_abc_d;
				best_case = 1;
				top_score = e_abc_d;
			}
			if (best_score + ERROR_TOLERANCE < cde_a_b + ce_ab_d + e_abd_c) {
				best_score = cde_a_b + ce_ab_d + e_abd_c;
				best_case = 2;
				top_score = e_abd_c;
			}
			//E|AB|CD
			switchSubtreeGet(c, 0, 2);
		}
		if (leafId(heavy(v)) == -1){
			int a = heavy(heavy(v)), b = light(heavy(v));
			//BE|A|CD
			switchSubtreeSet(b, 1, 0);
			tripScoreSet();
			//E|A|BCD
			switchSubtreeSet(b, 0, 2);
			tripScoreSet();
			//E|ACD|B
			switchSubtreeSet(light(v), 2, 1);
			tripScoreSet();
			//AE|B|CD
			switchSubtreeSet(light(v), 1, 2);
			switchSubtreeSet(a, 1, 0);
			switchSubtreeSet(b, 2, 1);
			tripScoreSet();
			//E|AB|CD
			switchSubtreeSet(a, 0, 1);
			
			//BE|A|CD
			switchSubtreeGet(b, 1, 0);
			score_t be_a_cd = tripScoreGet();
			//E|A|BCD
			switchSubtreeGet(b, 0, 2);
			score_t e_a_bcd = tripScoreGet();
			//E|ACD|B
			switchSubtreeGet(light(v), 2, 1);
			score_t e_acd_b = tripScoreGet();
			//AE|B|CD
			switchSubtreeGet(light(v), 1, 2);
			switchSubtreeGet(a, 1, 0);
			switchSubtreeGet(b, 2, 1);
			score_t ae_b_cd = tripScoreGet();
			if (best_score + ERROR_TOLERANCE < abe_c_d + ae_b_cd + e_a_bcd) {
				best_score = abe_c_d + ae_b_cd + e_a_bcd;
				best_case = 3;
				top_score = e_a_bcd;
			}
			if (best_score + ERROR_TOLERANCE < abe_c_d + be_a_cd + e_acd_b) {
				best_score = abe_c_d + be_a_cd + e_acd_b;
				best_case = 4;
				top_score = e_acd_b;
			}
			//E|AB|CD
			switchSubtreeGet(a, 0, 1);
		}
		//E|ABCD|-
		switchSubtree(light(v), 2, 1);
		return nnPerform(v, best_case, top_score);
	}
	
	score_t nnMove(int v){
		if (leafId(v) != -1){
			tripUpdateSet(1, leafId(v));
			tripUpdateGet(1, leafId(v));
			return 0;
		}
		else{
			//ABCDE|-|-
			score_t abe_c_d = nnMove(light(v));
			//ABE|CD|-
			
			score_t ae_b_cd = 0;
			if (leafId(heavy(v)) == -1){
				int b = light(heavy(v));
				switchSubtreeSet(b, 0, 2);
				//AE|CD|B
				tripScoreSet();
				switchSubtreeSet(b, 2, 0);
				//ABE|CD|-
				
				switchSubtreeGet(b, 0, 2);
				//AE|CD|B
				ae_b_cd = tripScoreGet();
				switchSubtreeGet(b, 2, 0);
				//ABE|CD|-
			}

			switchSubtree(light(v), 1, 0);
			//ABCDE|-|-
			score_t cde_a_b = nnMove(heavy(v));
			//CDE|AB|-
			switchSubtree(light(v), 0, 2);
			//E|AB|CD
			score_t e_ab_cd = tripScore();
			int best_case = 0;
			score_t best_score = abe_c_d + cde_a_b + e_ab_cd, top_score = e_ab_cd;
			
			if (leafId(light(v)) == -1){
				int c = heavy(light(v)), d = light(light(v));
				switchSubtreeSet(d, 2, 0);
				//DE|AB|C
				tripScoreSet();
				switchSubtreeSet(d, 0, 1);
				//E|ABD|C
				tripScoreSet();
				switchSubtreeSet(d, 1, 2);
				//E|AB|CD
				switchSubtreeSet(c, 2, 1);
				//E|ABC|D
				tripScoreSet();
				switchSubtreeSet(c, 1, 0);
				//CE|AB|D
				tripScoreSet();
				switchSubtreeSet(c, 0, 2);
				//E|AB|CD

				switchSubtreeGet(d, 2, 0);
				//DE|AB|C
				score_t de_ab_c = tripScoreGet();
				switchSubtreeGet(d, 0, 1);
				//E|ABD|C
				score_t e_abd_c = tripScoreGet();
				switchSubtreeGet(d, 1, 2);
				//E|AB|CD
				switchSubtreeGet(c, 2, 1);
				//E|ABC|D
				score_t e_abc_d = tripScoreGet();
				switchSubtreeGet(c, 1, 0);
				//CE|AB|D
				score_t ce_ab_d = tripScoreGet();
				
				if (best_score + ERROR_TOLERANCE < cde_a_b + de_ab_c + e_abc_d) {
					best_score = cde_a_b + de_ab_c + e_abc_d;
					best_case = 1;
					top_score = e_abc_d;
				}
				if (best_score + ERROR_TOLERANCE < cde_a_b + ce_ab_d + e_abd_c) {
					best_score = cde_a_b + ce_ab_d + e_abd_c;
					best_case = 2;
					top_score = e_abd_c;
				}
				switchSubtreeGet(c, 0, 2);
				//E|AB|CD
			}
			if (leafId(heavy(v)) == -1){
				int b = light(heavy(v));
				switchSubtreeSet(b, 1, 2);
				//E|A|BCD
				tripScoreSet();
				switchSubtreeSet(b, 2, 0);
				//BE|A|CD
				tripScoreSet();
				switchSubtreeSet(b, 0, 2);
				//E|A|BCD
				switchSubtreeSet(light(v), 2, 1);
				//E|ACD|B
				tripScoreSet();
				switchSubtreeSet(b, 2, 1);
				//E|ABCD|-
				switchSubtreeSet(light(v), 1, 2);
				//E|AB|CD

				switchSubtreeGet(b, 1, 2);
				//E|A|BCD
				score_t e_a_bcd = tripScoreGet();
				switchSubtreeGet(b, 2, 0);
				//BE|A|CD
				score_t be_a_cd = tripScoreGet();
				switchSubtreeGet(b, 0, 2);
				//E|A|BCD
				switchSubtreeGet(light(v), 2, 1);
				//E|ACD|B
				score_t e_acd_b = tripScoreGet();
				if (best_score + ERROR_TOLERANCE < abe_c_d + ae_b_cd + e_a_bcd) {
					best_score = abe_c_d + ae_b_cd + e_a_bcd;
					best_case = 3;
					top_score = e_a_bcd;
				}
				if (best_score + ERROR_TOLERANCE < abe_c_d + be_a_cd + e_acd_b) {
					best_score = abe_c_d + be_a_cd + e_acd_b;
					best_case = 4;
					top_score = e_acd_b;
				}
				switchSubtreeGet(b, 2, 1);
				//E|ABCD|-
				switchSubtreeGet(light(v), 1, 2);
				//E|AB|CD
			}
			switchSubtree(light(v), 2, 1);
			//E|ABCD|-
			return nnPerform(v, best_case, top_score);
		}
	}
	
	bool nnLocal2(int v, score_t& e_ab_cd) {
		//E|AB|CD
		if (leafId(v) != -1) return false;
		if (rNN <= 0) return false;
		rNN--;

		if (leafId(heavy(v)) == -1) {
			switchSubtreeSet(light(v), 2, 0);
			//CDE|AB|-
			switchSubtreeSet(light(heavy(v)), 1, 2);
			//CDE|A|B
			tripScoreSet();
			switchSubtreeSet(light(heavy(v)), 2, 1);
			//CDE|AB|-
			switchSubtreeSet(light(v), 0, 2);
			//E|AB|CD

			switchSubtreeGet(light(v), 2, 0);
			//CDE|AB|-
			switchSubtreeGet(light(heavy(v)), 1, 2);
			//CDE|A|B
			score_t cde_a_b = tripScoreGet();
			switchSubtreeGet(light(heavy(v)), 2, 1);
			//CDE|AB|-
			switchSubtreeGet(light(v), 0, 2);
			//E|AB|CD

			if (nnCheckHeavy(v, e_ab_cd, cde_a_b)) {
				nnLocal2(v, e_ab_cd);
				return true;
			}
		}

		if (leafId(light(v)) == -1) {
			switchSubtreeSet(heavy(v), 1, 0);
			//ABE|-|CD
			switchSubtreeSet(light(light(v)), 2, 1);
			//ABE|D|C
			tripScoreSet();
			switchSubtreeSet(light(light(v)), 1, 2);
			//ABE|-|CD
			switchSubtreeSet(heavy(v), 0, 1);
			//E|AB|CD

			switchSubtreeGet(heavy(v), 1, 0);
			//ABE|-|CD
			switchSubtreeGet(light(light(v)), 2, 1);
			//ABE|D|C
			score_t abe_c_d = tripScoreGet();
			switchSubtreeGet(light(light(v)), 1, 2);
			//ABE|-|CD
			switchSubtreeGet(heavy(v), 0, 1);
			//E|AB|CD

			if (nnCheckLight(v, e_ab_cd, abe_c_d)) {
				nnLocal2(v, e_ab_cd);
				return true;
			}
		}
		
		return false;
	}

	bool nnCheckLight(int v, score_t& e_ab_cd, score_t& abe_c_d) {
		//E|AB|CD
		if (leafId(light(v)) == -1) {
			int u = light(v), c = heavy(u), d = light(u), ab = heavy(v);
			// int c = heavy(light(v)), d = light(light(v));
			switchSubtreeSet(d, 2, 0);
			//DE|AB|C
			tripScoreSet();
			switchSubtreeSet(d, 0, 1);
			//E|ABD|C
			tripScoreSet();
			switchSubtreeSet(d, 1, 2);
			//E|AB|CD
			switchSubtreeSet(c, 2, 1);
			//E|ABC|D
			tripScoreSet();
			switchSubtreeSet(c, 1, 0);
			//CE|AB|D
			tripScoreSet();
			switchSubtreeSet(c, 0, 2);
			//E|AB|CD

			switchSubtreeGet(d, 2, 0);
			//DE|AB|C
			score_t de_ab_c = tripScoreGet();
			switchSubtreeGet(d, 0, 1);
			//E|ABD|C
			score_t e_abd_c = tripScoreGet();
			switchSubtreeGet(d, 1, 2);
			//E|AB|CD
			switchSubtreeGet(c, 2, 1);
			//E|ABC|D
			score_t e_abc_d = tripScoreGet();
			switchSubtreeGet(c, 1, 0);
			//CE|AB|D
			score_t ce_ab_d = tripScoreGet();
			switchSubtreeGet(c, 0, 2);
			//E|AB|CD

			score_t best_score = e_ab_cd + abe_c_d;
			int best_case = 0;
			if (best_score + ERROR_TOLERANCE < de_ab_c + e_abc_d) {
				best_score = de_ab_c + e_abc_d;
				best_case = 1;
			}
			if (best_score + ERROR_TOLERANCE < ce_ab_d + e_abd_c) {
				best_score = ce_ab_d + e_abd_c;
				best_case = 2;
			}

			if (best_case == 1) {
				parent(ab) = u; heavy(u) = ab;
				parent(c) = u; light(u) = c;
				parent(u) = v; heavy(v) = u;
				parent(d) = v; light(v) = d;
				leafCnt(u) = leafCnt(heavy(u)) + leafCnt(light(u));
				e_ab_cd = e_abc_d;
				abe_c_d = de_ab_c;

				switchSubtree(d, 2, 0);
				//DE|AB|C
				nnLocal2(u, abe_c_d);

				switchSubtreeSet(d, 0, 2);
				//E|AB|CD
				switchSubtreeSet(c, 2, 1);
				//E|ABC|D

				switchSubtreeGet(d, 0, 2);
				//E|AB|CD
				switchSubtreeGet(c, 2, 1);
				//E|ABC|D
				return true;
			}
			if (best_case == 2) {
				switchSubtree(d, 2, 1);
				//E|ABD|C
				parent(ab) = u; heavy(u) = ab;
				parent(d) = u; light(u) = d;
				parent(u) = v; heavy(v) = u;
				parent(c) = v; light(v) = c;
				leafCnt(u) = leafCnt(heavy(u)) + leafCnt(light(u));
				e_ab_cd = e_abd_c;
				abe_c_d = ce_ab_d;

				switchSubtree(c, 2, 0);
				//CE|AB|D
				nnLocal2(u, abe_c_d);

				switchSubtreeSet(c, 0, 2);
				//E|AB|CD
				switchSubtreeSet(d, 2, 1);
				//E|ABD|C

				switchSubtreeGet(c, 0, 2);
				//E|AB|CD
				switchSubtreeGet(d, 2, 1);
				//E|ABD|C
				return true;
			}
		}
		return false;
	}

	score_t nnPrecheckHeavy(int v) {
		//ABE|CD|-
		score_t ae_b_cd = 0;
		if (leafId(heavy(v)) == -1) {
			int b = light(heavy(v));
			switchSubtreeSet(b, 0, 2);
			//AE|CD|B
			tripScoreSet();
			switchSubtreeSet(b, 2, 0);
			//ABE|CD|-

			switchSubtreeGet(b, 0, 2);
			//AE|CD|B
			ae_b_cd = tripScoreGet();
			switchSubtreeGet(b, 2, 0);
			//ABE|CD|-
		}

		return ae_b_cd;
	}
	
	bool nnCheckHeavy(int v, score_t& e_ab_cd, score_t& cde_a_b, bool precheck = false, score_t ae_b_cd = 0) {
		//E|AB|CD
		if (leafId(heavy(v)) == -1) {
			int u = heavy(v), a = heavy(u), b = light(u), cd = light(v);
			if (!precheck) {
				switchSubtreeSet(a, 1, 0);
				//AE|B|CD
				tripScoreSet();
				switchSubtreeSet(a, 0, 1);
				//E|AB|CD

				switchSubtreeGet(a, 1, 0);
				//AE|B|CD
				ae_b_cd = tripScoreGet();
				switchSubtreeGet(a, 0, 1);
				//E|AB|CD
			}

			switchSubtreeSet(b, 1, 2);
			//E|A|BCD
			tripScoreSet();
			switchSubtreeSet(b, 2, 0);
			//BE|A|CD
			tripScoreSet();
			switchSubtreeSet(b, 0, 2);
			//E|A|BCD
			switchSubtreeSet(cd, 2, 1);
			//E|ACD|B
			tripScoreSet();
			switchSubtreeSet(b, 2, 1);
			//E|ABCD|-
			switchSubtreeSet(cd, 1, 2);
			//E|AB|CD

			switchSubtreeGet(b, 1, 2);
			//E|A|BCD
			score_t e_a_bcd = tripScoreGet();
			switchSubtreeGet(b, 2, 0);
			//BE|A|CD
			score_t be_a_cd = tripScoreGet();
			switchSubtreeGet(b, 0, 2);
			//E|A|BCD
			switchSubtreeGet(cd, 2, 1);
			//E|ACD|B
			score_t e_acd_b = tripScoreGet();
			switchSubtreeGet(b, 2, 1);
			//E|ABCD|-
			switchSubtreeGet(cd, 1, 2);
			//E|AB|CD

			score_t best_score = e_ab_cd + cde_a_b;
			int best_case = 0;
			if (best_score + ERROR_TOLERANCE < ae_b_cd + e_a_bcd) {
				best_score = ae_b_cd + e_a_bcd;
				best_case = 3;
			}
			if (best_score + ERROR_TOLERANCE < be_a_cd + e_acd_b) {
				best_score = be_a_cd + e_acd_b;
				best_case = 4;
			}

			if (best_case == 3) {
				parent(b) = u; heavy(u) = b;
				parent(cd) = u; light(u) = cd;
				parent(a) = v; heavy(v) = a;
				parent(u) = v; light(v) = u;
				leafCnt(u) = leafCnt(heavy(u)) + leafCnt(light(u));
				e_ab_cd = e_a_bcd;
				cde_a_b = ae_b_cd;

				switchSubtree(a, 1, 0);
				//AE|B|CD
				nnLocal2(u, cde_a_b);

				switchSubtreeSet(a, 0, 1);
				//E|AB|CD
				switchSubtreeSet(b, 1, 2);
				//E|A|BCD

				switchSubtreeGet(a, 0, 1);
				//E|AB|CD
				switchSubtreeGet(b, 1, 2);
				//E|A|BCD
				return true;
			}
			if (best_case == 4) {
				parent(a) = u; heavy(u) = a;
				parent(cd) = u; light(u) = cd;
				parent(u) = v; heavy(v) = u;
				parent(b) = v; light(v) = b;
				leafCnt(u) = leafCnt(heavy(u)) + leafCnt(light(u));
				e_ab_cd = e_acd_b;
				cde_a_b = be_a_cd;

				switchSubtree(b, 1, 0);
				//BE|A|CD
				nnLocal2(u, cde_a_b);
				
				switchSubtreeSet(cd, 2, 1);
				//BE|ACD|-
				switchSubtreeSet(b, 0, 2);
				//E|ACD|B

				switchSubtreeGet(cd, 2, 1);
				//BE|ACD|-
				switchSubtreeGet(b, 0, 2);
				//E|ACD|B

				return true;
			}
		}
		return false;
	}

	score_t nnMove2(int v) {
		if (leafId(v) != -1) {
			tripUpdateSet(1, leafId(v));
			tripUpdateGet(1, leafId(v));
			return 0;
		}
		else {
			//ABCDE|-|-
			score_t abe_c_d = nnMove2(light(v));
			//ABE|CD|-

			score_t ae_b_cd = nnPrecheckHeavy(v);

			switchSubtree(light(v), 1, 0);
			//ABCDE|-|-
			score_t cde_a_b = nnMove2(heavy(v));
			//CDE|AB|-
			switchSubtree(light(v), 0, 2);
			//E|AB|CD
			score_t e_ab_cd = tripScore();

			if (nnCheckHeavy(v, e_ab_cd, cde_a_b, true, ae_b_cd)) nnLocal2(v, e_ab_cd);
			else if (nnCheckLight(v, e_ab_cd, abe_c_d)) nnLocal2(v, e_ab_cd);

			switchSubtree(light(v), 2, 1);
			//E|ABCD|-

			return e_ab_cd;
		}
	}

	void nnMove(){
		if (ROUND_NN >= 0) {
			tripUpdateSet(0, rootLeafId);
			switchSubtreeSet(rootNodeId, -1, 0);
			tripUpdateGet(0, rootLeafId);
			switchSubtreeGet(rootNodeId, -1, 0);
			nnMove2(rootNodeId);
			LOG << "#NNI moves:" << ROUND_NN - rNN << "/" << ROUND_NN << endl;
		}
	}
	
	void tripHashGeneratorSet(int v){
		if (leafId(v) != -1){
			tripUpdateSet(1, leafId(v));
		}
		else{
			//AC|B|-
			tripHashGeneratorSet(light(v));
			//ABC|-|-
			switchSubtreeSet(light(v), 1, 0);
			//BC|A|-
			tripHashGeneratorSet(heavy(v));
			//C|A|B
			switchSubtreeSet(light(v), 0, 2);
			tripScoreSet();
			//C|AB|-
			switchSubtreeSet(light(v), 2, 1);
		}
		
	}

	pair<bool, hash_t> tripHashGeneratorGet(int v){
		if (leafId(v) != -1){
			tripUpdateGet(1, leafId(v));
			return {leafId(v) == 0, taxonHash[leafId(v)]};
		}
		else{
			//AC|B|-
			pair<bool, hash_t> hLight = tripHashGeneratorGet(light(v));
			//ABC|-|-
			switchSubtreeGet(light(v), 1, 0);
			//BC|A|-
			pair<bool, hash_t> hHeavy = tripHashGeneratorGet(heavy(v));
			//C|A|B
			switchSubtreeGet(light(v), 0, 2);
			score_t s = tripScoreGet();
			//C|AB|-
			switchSubtreeGet(light(v), 2, 1);
			if (hHeavy.first) tripHash.push_back(make_tuple(hLight.second, -hHeavy.second-hLight.second, s));
			else if (hLight.first) tripHash.push_back(make_tuple(hHeavy.second, -hHeavy.second-hLight.second, s));
			else tripHash.push_back(make_tuple(hHeavy.second, hLight.second, s));
			return {hHeavy.first || hLight.first, hHeavy.second + hLight.second};
		}
		
	}
	
	void run(){
		defaultInitializer();
		for (; orderId < order.size(); orderId++){
			place(order[orderId]);
		}
		nnMove();
		tripUpdateSet(0, rootLeafId);
		switchSubtreeSet(rootNodeId, -1, 0);
		tripHashGeneratorSet(rootNodeId);
		tripUpdateGet(0, rootLeafId);
		switchSubtreeGet(rootNodeId, -1, 0);
		tripHashGeneratorGet(rootNodeId);
	}

	string printSubtree(int v){
		if (leafId(v) != -1) return to_string(leafId(v));
		return string("(") + printSubtree(heavy(v)) + "," + printSubtree(light(v)) + ")";
	}
	
	string printTree(){
		return string("(") + printSubtree(rootNodeId) + "," + to_string(rootLeafId) + ");";
	}
	
	string printSubtree(int v, const vector<string> &names){
		if (leafId(v) != -1) return names[leafId(v)];
		return string("(") + printSubtree(heavy(v), names) + "," + printSubtree(light(v), names) + ")";
	}
	
	string printTree(const vector<string> &names){
		return string("(") + printSubtree(rootNodeId, names) + "," + names[rootLeafId] + ");";
	}
};

struct ConstrainedOptimizationAlgorithm{
	struct Node{
		Node(int leafId, hash_t hash): hash(hash), leafId(leafId) {}
		
		hash_t hash, bestChild = 0;
		score_t bestScore = 0;
		int roundId = 0, leafId;
		unordered_map<hash_t, tuple<int, int, score_t>, Hasher> children;
	};
	
	const int ntaxa;
	TripartitionInitializer& tripInit;
	const vector<string> names;
	vector<Node> nodes;
	unordered_map<hash_t, int, Hasher> hash;
	vector<hash_t> taxonHash;
	uniform_int_distribution<hash_t> randomHash;
	uniform_real_distribution<double> randP;
	int roundId = 0;
	ThreadPool &TP;
	int ROUND_NN = -1;

	shared_ptr<AnnotatedTree> annotTree;
	double tempSibling, tempRoot;

	ConstrainedOptimizationAlgorithm(const int ntaxa, TripartitionInitializer &tripInit, const vector<string> &names, ThreadPool &TP, int ROUND_NN, const int seed = rand()):
			ntaxa(ntaxa), tripInit(tripInit), names(names), TP(TP), ROUND_NN(ROUND_NN) {
		taxonHash.push_back(0);
		for (int i = 1; i < ntaxa; i++){
			hash_t r = randomHash(RND_GENERATOR);
			taxonHash.push_back(r);
			hash[r] = nodes.size();
			nodes.emplace_back(i, r);
			taxonHash[0] -= r;
		}
	}
	
	ConstrainedOptimizationAlgorithm(const ConstrainedOptimizationAlgorithm &alg): 
			ntaxa(alg.ntaxa), tripInit(alg.tripInit), names(alg.names), nodes(alg.nodes), hash(alg.hash), taxonHash(alg.taxonHash), TP(alg.TP), ROUND_NN(alg.ROUND_NN) {}
	
	int subsampleSubtree(int v, PlacementAlgorithm &pAlg, const unordered_set<int> &selected){
		if (nodes[v].leafId != -1){
			if (selected.count(nodes[v].leafId)){
				if (pAlg.rootLeafId == -1){
					pAlg.rootLeafId = nodes[v].leafId;
					return -1;
				}
				else {
					pAlg.nodes.emplace_back(nodes[v].leafId);
					return pAlg.nodes.size() - 1;
				}
			}
			else {
				pAlg.order.push_back(nodes[v].leafId);
				return -1;
			}
		}
		else {
			tuple<int, int, score_t> t = nodes[v].children[nodes[v].bestChild];
			if (pAlg.rootLeafId != -1){
				int nodeL = subsampleSubtree(get<0>(t), pAlg, selected);
				int nodeR = subsampleSubtree(get<1>(t), pAlg, selected);
				if (nodeL != -1){
					if (nodeR != -1) return pAlg.makeNode(nodeL, nodeR, false);
					else return nodeL;
				}
				else return nodeR;
			}
			else {
				int nodeL = subsampleSubtree(get<0>(t), pAlg, selected);
				if (pAlg.rootLeafId == -1) return subsampleSubtree(get<1>(t), pAlg, selected);
				else if (nodeL == -1) return pAlg.rootNodeId = subsampleSubtree(get<1>(t), pAlg, selected);
				else {
					int nodeR = subsampleSubtree(get<1>(t), pAlg, selected);
					if (nodeR == -1) return nodeL;
					else {
						pAlg.makeNode(nodeL, nodeR, true);
						return nodeR;
					}
				}
			}
		}
	}
	
	void createPlacementAlgorithm(PlacementAlgorithm &pAlg, const unordered_set<int> &selected){
		if (roundId != 0){
			if (selected.count(0)){
				pAlg.rootLeafId = 0;
				pAlg.rootNodeId = subsampleSubtree(hash[-taxonHash[0]], pAlg, selected);
			}
			else {
				pAlg.order.push_back(0);
				subsampleSubtree(hash[-taxonHash[0]], pAlg, selected);
			}
		}
		else{
			for (int i = 0; i < ntaxa; i++){
				pAlg.order.push_back(i);
			}
		}
		shuffle(pAlg.order.begin(), pAlg.order.end(), RND_GENERATOR);
	}
	
	int subsampleSubtree(int v, PlacementAlgorithm &pAlg, double subsampleRate){
		if (nodes[v].leafId != -1){
			if (randP(RND_GENERATOR) < subsampleRate){
				if (pAlg.rootLeafId == -1){
					pAlg.rootLeafId = nodes[v].leafId;
					return -1;
				}
				else {
					pAlg.nodes.emplace_back(nodes[v].leafId);
					return pAlg.nodes.size() - 1;
				}
			}
			else {
				pAlg.order.push_back(nodes[v].leafId);
				return -1;
			}
		}
		else {
			tuple<int, int, score_t> t = nodes[v].children[nodes[v].bestChild];
			if (pAlg.rootLeafId != -1){
				int nodeL = subsampleSubtree(get<0>(t), pAlg, subsampleRate);
				int nodeR = subsampleSubtree(get<1>(t), pAlg, subsampleRate);
				if (nodeL != -1){
					if (nodeR != -1) return pAlg.makeNode(nodeL, nodeR, false);
					else return nodeL;
				}
				else return nodeR;
			}
			else {
				int nodeL = subsampleSubtree(get<0>(t), pAlg, subsampleRate);
				if (pAlg.rootLeafId == -1) return subsampleSubtree(get<1>(t), pAlg, subsampleRate);
				else if (nodeL == -1) return pAlg.rootNodeId = subsampleSubtree(get<1>(t), pAlg, subsampleRate);
				else {
					int nodeR = subsampleSubtree(get<1>(t), pAlg, subsampleRate);
					if (nodeR == -1) return nodeL;
					else {
						pAlg.makeNode(nodeL, nodeR, true);
						return nodeR;
					}
				}
			}
		}
	}
	
	void createPlacementAlgorithm(PlacementAlgorithm &pAlg, double subsampleRate){
		if (roundId != 0){
			if (randP(RND_GENERATOR) < subsampleRate){
				pAlg.rootLeafId = 0;
				//pAlg.trip.addTotal(0);
				pAlg.rootNodeId = subsampleSubtree(hash[-taxonHash[0]], pAlg, subsampleRate);
			}
			else {
				pAlg.order.push_back(0);
				subsampleSubtree(hash[-taxonHash[0]], pAlg, subsampleRate);
			}
		}
		else{
			for (int i = 0; i < ntaxa; i++){
				pAlg.order.push_back(i);
			}
		}
		shuffle(pAlg.order.begin(), pAlg.order.end(), RND_GENERATOR);
	}
	
	int guideSubtree(PlacementAlgorithm &pAlg, const string &tree, const unordered_map<string, int> &name2id, int &i, unordered_set<int> &added){
		int ret;
		if (tree[i] != '('){
			string s;
			for (; tree[i] != ':' && tree[i] != ')' && tree[i] != ',' && tree[i] != ';'; i++){
				if (tree[i] != '\"' && tree[i] != '\'') s += tree[i];
			}
			if (name2id.count(s) == 0){
				cerr << "Unknow species '" << s << "' in tree '" << tree << "' not found in the input! \n";
				exit(0);
			}
			int id = name2id.at(s);
			//pAlg.trip.addTotal(id);
			added.insert(id);
			if (pAlg.rootLeafId == -1){
				pAlg.rootLeafId = id;
				ret = -1;
			}
			else {
				pAlg.nodes.emplace_back(id);
				ret = pAlg.nodes.size() - 1;
			}
		}
		else {
			i++;
			if (pAlg.rootLeafId != -1){
				int nodeL = guideSubtree(pAlg, tree, name2id, i, added);
				int nodeR = guideSubtree(pAlg, tree, name2id, i, added);
				if (nodeL != -1){
					if (nodeR != -1) ret = pAlg.makeNode(nodeL, nodeR, false);
					else ret = nodeL;
				}
				else ret = nodeR;
			}
			else {
				int nodeL = guideSubtree(pAlg, tree, name2id, i, added);
				if (pAlg.rootLeafId == -1) ret = guideSubtree(pAlg, tree, name2id, i, added);
				else if (nodeL == -1) ret = pAlg.rootNodeId = guideSubtree(pAlg, tree, name2id, i, added);
				else {
					int nodeR = guideSubtree(pAlg, tree, name2id, i, added);
					if (nodeR == -1) ret = nodeL;
					else {
						pAlg.makeNode(nodeL, nodeR, true);
						ret = nodeR;
					}
				}
			}
			//while (tree[i] != ')' && tree[i] != ',' && tree[i] != ';') i++;
			if (tree[i-1] == ',') {
				LOG << "Warning: Polytomy detected! Resolving arbitrarily!\n";
				//i++;
				int nodeL = ret;
				if (pAlg.rootLeafId == -1) ret = guideSubtree(pAlg, tree, name2id, i, added);
				else if (nodeL == -1) ret = pAlg.rootNodeId = guideSubtree(pAlg, tree, name2id, i, added);
				else {
					int nodeR = guideSubtree(pAlg, tree, name2id, i, added);
					if (nodeR == -1) ret = nodeL;
					else {
						pAlg.makeNode(nodeL, nodeR, true);
						ret = nodeR;
					}
				}
			}
		}
		while (tree[i] != ')' && tree[i] != ',' && tree[i] != ';') i++;
		i++;
		return ret;
	}

	void addGuideTree(const string &tree, const unordered_map<string, int> &name2id, int rnn){
		unordered_set<int> added;
		PlacementAlgorithm pAlg(taxonHash, tripInit, TP, rnn);
		int i = 0;
		guideSubtree(pAlg, tree, name2id, i, added);
		for (int i = 0; i < ntaxa; i++){
			if (added.count(i) == 0) pAlg.order.push_back(i);
		}
		shuffle(pAlg.order.begin(), pAlg.order.end(), RND_GENERATOR);
		pAlg.run();
		LOG << pAlg.printTree(names) << endl;
		addTripartitions(pAlg.tripHash);
	}
	
	int hash2node(hash_t h){
		if (hash.count(h) == 0){
			hash[h] = nodes.size();
			nodes.emplace_back(-1, h);
		}
		return hash[h];
	}
	
	void addTripartitions(const vector<tuple<hash_t, hash_t, score_t> > tripHash){
		for (const tuple<hash_t, hash_t, score_t> &t: tripHash){
			int left = hash2node(get<0>(t)), right = hash2node(get<1>(t)), parent = hash2node(get<0>(t) + get<1>(t));
			nodes[parent].children[get<0>(t)] = make_tuple(left, right, get<2>(t));
			nodes[parent].children[get<1>(t)] = make_tuple(left, right, get<2>(t));
		}
	}
	
	void optimalTreeDP(int v){
		if (nodes[v].roundId == roundId) return;
		for (pair<hash_t, tuple<int, int, score_t> > child: nodes[v].children){
			optimalTreeDP(get<0>(child.second));
			optimalTreeDP(get<1>(child.second));
			if (nodes[v].roundId < roundId || nodes[v].bestScore < nodes[get<0>(child.second)].bestScore + nodes[get<1>(child.second)].bestScore + get<2>(child.second)){
				nodes[v].roundId = roundId;
				nodes[v].bestScore = nodes[get<0>(child.second)].bestScore + nodes[get<1>(child.second)].bestScore + get<2>(child.second);
				nodes[v].bestChild = child.first;
			}
		}
		nodes[v].roundId = roundId;
	}
	
	score_t computeOptimalTree(){
		roundId++;
		optimalTreeDP(hash2node(-taxonHash[0]));
		return nodes[hash[-taxonHash[0]]].bestScore;
	}
	
	string printOptimalSubtree(int v){
		if (nodes[v].leafId != -1) return names[nodes[v].leafId];
		tuple<int, int, score_t> t = nodes[v].children[nodes[v].bestChild];
		return string("(") + printOptimalSubtree(get<0>(t)) + "," + printOptimalSubtree(get<1>(t)) + ")";
	}
	
	string printOptimalTree(){
		return string("(") + printOptimalSubtree(hash[-taxonHash[0]]) + "," + names[0] + ");";
	}
	
	string printOptimalSubtreeWithScore(int v){
		if (nodes[v].leafId != -1) return names[nodes[v].leafId];
		tuple<int, int, score_t> t = nodes[v].children[nodes[v].bestChild];
		return string("(") + printOptimalSubtreeWithScore(get<0>(t)) + "," + printOptimalSubtreeWithScore(get<1>(t)) + "):" + to_string(get<2>(t));
	}
	
	string printOptimalTreeWithScore(){
		return string("(") + printOptimalSubtreeWithScore(hash[-taxonHash[0]]) + "," + names[0] + ");";
	}
	
	pair<score_t, string> run(int nJobs = 1, double subsampleRate = 0, bool allowTwoStepRun = true){
		if (nJobs == 0) return {computeOptimalTree(), printOptimalTree()};
		if (allowTwoStepRun && ntaxa >= 100) return twoStepRun(nJobs);
		for (int i = 0; i < nJobs; i++) {
			PlacementAlgorithm pAlg(taxonHash, tripInit, TP, ROUND_NN);
			createPlacementAlgorithm(pAlg, subsampleRate);
			pAlg.run();
			LOG << pAlg.printTree(names) << endl;
			addTripartitions(pAlg.tripHash);
		}
		return {computeOptimalTree(), printOptimalTree()};
	}

	pair<score_t, string> constrainedRun(int nJobs, const string tree, const unordered_map<string, int> &name2id){
		for (int i = 0; i < nJobs; i++) {
			addGuideTree(tree, name2id, -1);
		}
		return {computeOptimalTree(), printOptimalTree()};
	}
	
	static bool placeNode(PlacementAlgorithm &pAlg, vector<PlacementAlgorithm::Node> &pNodes, int &rootNodeId, int i, unordered_map<int, int> &nodeRemap){
		int v = pAlg.locateBranchAdding(i);
		if (nodeRemap.count(v)){
			pAlg.addSibling(v, i);
			int w = pAlg.parent(v), u = (pAlg.light(w) == v) ? pAlg.heavy(w) : pAlg.light(w);
			int mV = nodeRemap[v], mW = pNodes.size(), mU = pNodes.size() + 1;
			pNodes.emplace_back(-1);
			pNodes.emplace_back(i);
			nodeRemap[w] = mW;
			nodeRemap[u] = mU;
			if (pAlg.rootNodeId == w) rootNodeId = mW;
			else {
				if (pNodes[pNodes[mV].parent].heavy == mV) pNodes[pNodes[mV].parent].heavy = mW;
				else pNodes[pNodes[mV].parent].light = mW;
			}
			pNodes[mW].parent = pNodes[mV].parent;
			pNodes[mV].parent = mW;
			pNodes[mU].parent = mW;
			pNodes[mW].heavy = mV;
			pNodes[mW].light = mU;
			return true;
		}
		else {
			pAlg.tripUpdateSet(-1, i);
			pAlg.tripUpdateGet(-1, i);
			return false;
		}
	}
	
	void twoStepWorkflow(PlacementAlgorithm &pAlg){
		int N = ntaxa, rId = roundId;
		int n = sqrt(N) * log2(names.size()) / 4;
		vector<int> order;
		if (rId == 0){
			ConstrainedOptimizationAlgorithm alg(*this);
			for (int i = 1; i < N; i++) order.push_back(i);
			shuffle(order.begin(), order.end(), RND_GENERATOR);
			order.push_back(order[0]);
			order[0] = 0;
			hash_t hashsum = 0;
			for (int i = 1; i < n; i++) hashsum += taxonHash[order[i]];
			
			pAlg.taxonHash[0] = -hashsum;
			alg.taxonHash[0] = -hashsum;
			for (int r = 0; r < n; r++){
				LOG << "Guide Tree " << r << "/" << n << endl;
				for (int i = 0; i < n; i++) pAlg.order.push_back(order[i]);
				shuffle(pAlg.order.begin(), pAlg.order.end(), RND_GENERATOR);
				pAlg.run();
				alg.addTripartitions(pAlg.tripHash);
				for (int i = 0; i < n; i++) pAlg.tripUpdateSet(-1, order[i]);
				for (int i = 0; i < n; i++) pAlg.tripUpdateGet(-1, order[i]);
				pAlg.tripHash.clear();
				pAlg.order.clear();
				pAlg.nodes.clear();
				pAlg.rootNodeId = -1;
				pAlg.rootLeafId = -1;
				pAlg.orderId = 0;
				pAlg.rNN = ROUND_NN;
			}
			alg.computeOptimalTree();
			alg.createPlacementAlgorithm(pAlg, 1);
			pAlg.taxonHash[0] = taxonHash[0];
		}
		else{
			for (int i = 0; i < N; i++) order.push_back(i);
			shuffle(order.begin(), order.end(), RND_GENERATOR);
			unordered_set<int> selected;
			for (int i = 0; i < n; i++) selected.insert(order[i]);
			createPlacementAlgorithm(pAlg, selected);
			pAlg.nnMove();
		}
		{
			vector<int> abnormalOrder;
			vector<PlacementAlgorithm::Node> pNodes(pAlg.nodes);
			int rootNodeId = pAlg.rootNodeId;
			{
				vector<PlacementAlgorithm::Node> originalNodes(pAlg.nodes);
				int originalRootNodeId = pAlg.rootNodeId;
				
				vector<vector<int> > nodeBranch(pNodes.size());
				
				for (int i = n; i < N; i++) {
					nodeBranch[pAlg.locateBranch(order[i])].push_back(order[i]);
				}
				for (size_t b = 0; b < nodeBranch.size(); b++){
					pAlg.nodes = originalNodes;
					pAlg.rootNodeId = originalRootNodeId;
					vector<int> added;
					unordered_map<int, int> nodeRemap;
					nodeRemap[b] = b;
					for (int i: nodeBranch[b]) {
						if (placeNode(pAlg, pNodes, rootNodeId, i, nodeRemap)) added.push_back(i);
						else abnormalOrder.push_back(i);
					}
					for (int i: added) pAlg.tripUpdateSet(-1, i);
					for (int i: added) pAlg.tripUpdateGet(-1, i);
				}
			}
			pAlg.order = abnormalOrder;
			pAlg.nodes = pNodes;
			pAlg.rootNodeId = rootNodeId;
			pAlg.orderId = 0;
		}
		LOG << "Remaining: " << pAlg.order.size() << endl;
		pAlg.nnMove();
		pAlg.run();
	}
	
	pair<score_t, string> twoStepRun(int nJobs){
		LOG << "Use two-step algorithm!\n";
		
		for (int i = 0; i < nJobs; i++) {
			PlacementAlgorithm job(taxonHash, tripInit, TP, ROUND_NN);
			twoStepWorkflow(job);
			LOG << job.printTree(names) << endl;
			addTripartitions(job.tripHash);
		}
		
		return {computeOptimalTree(), printOptimalTree()};
	}

	inline static double gTest(double v1, double v2, double v3){
		if (v1 < max(v2, v3) - 1e-10) return 0;
		if (abs(v2 - v3) < 1e-10) {
			if (v1 - max(v2, v3) < 1e-10) return 0.5;
			else return 1;
		}
		double avg = (v1 + v2 + v3) / 3;
		double sd = sqrt((v1 * v1 + v2 * v2 + v3 * v3 - 3 * avg * avg) / 2);
		double G = (v1 - avg) / sd;
		double t = sqrt((3 * G * G) / (4 - 3 * G * G));
		double p = 0.5 - atan(t) / 3.14159265;
		return 1 - p * 3;
	}

#ifdef SOME_FORM_OF_SUPPORT
	void switchSubtree(Quadrupartition &quad, int v, int s, int t){
		if (nodes[v].leafId != -1) {
			quad.update(t, nodes[v].leafId);
		}
		else {
			tuple<int, int, score_t> c = nodes[v].children[nodes[v].bestChild];
			switchSubtree(quad, get<0>(c), s, t);
			switchSubtree(quad, get<1>(c), s, t);
		}	
	}

	/*
	string printOptimalSubtreeWithSupport(Quadrupartition &quad, int v, int u, int support, double lambda,
			unordered_map<int, tuple<array<double, 3>, array<double, 3>, string> > &qInfo, double weight, shared_ptr<AnnotatedTree::Node> node, bool toplevel = false){
		if (nodes[v].leafId != -1) {
			int i = nodes[v].leafId;
			array<double, 3> t;
			if (support == 3) qInfo[v] = make_tuple(t, t, names[i]);
			node->taxon = i;
			node->name = names[i];
			return names[i];
		}
		//r|u|c0c1|-
		string res;
		tuple<int, int, score_t> c = nodes[v].children[nodes[v].bestChild];
		//ru|c1|c0|-
		switchSubtree(quad, u, 1, 0);
		switchSubtree(quad, get<1>(c), 2, 1);
		string resL = printOptimalSubtreeWithSupport(quad, get<0>(c), get<1>(c), support, lambda, qInfo, weight, annotTree->addLeft(node));
		//ru|c0|c1|-
		switchSubtree(quad, get<0>(c), 2, 1);
		switchSubtree(quad, get<1>(c), 1, 2);
		string resR = printOptimalSubtreeWithSupport(quad, get<1>(c), get<0>(c), support, lambda, qInfo, weight, annotTree->addRight(node));
		//r|u|c1|c0
		switchSubtree(quad, u, 0, 1);
		switchSubtree(quad, get<0>(c), 1, 3);
		#ifdef CUSTOMIZED_ANNOTATION
		node->annot = quad.annotate();
		#endif
		array<double, 3> score = quad.score();
		//r|u|c0c1|-
		switchSubtree(quad, get<0>(c), 3, 2);
		score[0] /= weight; score[1] /= weight; score[2] /= weight;
		double tscore = score[0] + score[1] + score[2];
		#ifdef SUPPORT
		double support0, support1, support2;
		if (tscore < 1e-8){
			support0 = 1.0 / 3; support1 = 1.0 / 3; support2 = 1.0 / 3;
		}
		else if (tscore > 100 && score[0] - max(score[1], score[2]) > 5 * sqrt(tscore)){
			support0 = 1; support1 = 0; support2 = 0;
		}
		else if (tscore > 1000 && score[0] - max(score[1], score[2]) > 5 * sqrt((score[0] + max(score[1], score[2])) * 0.25)){
			support0 = 1; support1 = 0; support2 = 0;
		}
		else {
			double i0 = 1.0 - incbeta(score[0] + 1.0, tscore + lambda * 2 - score[0], 1.0 / 3.0);
			double i1 = 1.0 - incbeta(score[1] + 1.0, tscore + lambda * 2 - score[1], 1.0 / 3.0);
			double i2 = 1.0 - incbeta(score[2] + 1.0, tscore + lambda * 2 - score[2], 1.0 / 3.0);
			double lb0 = lgamma(score[0] + 1.0) + lgamma(tscore - score[0] + lambda * 2);
			double lb1 = lgamma(score[1] + 1.0) + lgamma(tscore - score[1] + lambda * 2);
			double lb2 = lgamma(score[2] + 1.0) + lgamma(tscore - score[2] + lambda * 2);
			support0 = i0 / (i0 + i1 * exp(log(2.0) * (score[1] - score[0]) + lb1 - lb0) + i2 * exp(log(2.0) * (score[2] - score[0]) + lb2 - lb0));
			support1 = i1 / (i1 + i0 * exp(log(2.0) * (score[0] - score[1]) + lb0 - lb1) + i2 * exp(log(2.0) * (score[2] - score[1]) + lb2 - lb1));
			support2 = i2 / (i2 + i1 * exp(log(2.0) * (score[1] - score[2]) + lb1 - lb2) + i0 * exp(log(2.0) * (score[0] - score[2]) + lb0 - lb2));
		}
		if (support == 1) res += to_string(support0);
		else {
			array<double, 3> p = {support0, support1, support2};
			res += "'[pp1=" + to_string(p[0]) + ";pp2=" + to_string(p[1]) + ";pp3=" + to_string(p[2]) + ";f1=" + to_string(score[0]) + ";f2=" + to_string(score[1]) + ";f3=" + to_string(score[2]) 
				+ ";q1=" + to_string(score[0] / tscore) + ";q2=" + to_string(score[1] / tscore) + ";q3=" + to_string(score[2] / tscore) + "]'";
			if (support == 3) qInfo[v] = make_tuple(score, p, get<2>(qInfo[get<0>(c)]) + "," + get<2>(qInfo[get<1>(c)]));
		}
		node->s = support0;
		#endif
		#ifdef LOCAL_BOOTSTRAP
		array<int, 3> bs = node->annot.bootstrap();
		double nbs = node->annot.bs.size() / 100.0;
		if (support == 1) res += formatBootstrap(bs[0] / nbs);
		else {
			res += "'[bs1=" + to_string(bs[0]) + ";bs2=" + to_string(bs[1]) + ";bs3=" + to_string(bs[2]) + 
				   ";s1=" + to_string(score[0]) + ";s2=" + to_string(score[1]) + ";s3=" + to_string(score[2]) +
				   ";q1=" + to_string(score[0] / tscore) + ";q2=" + to_string(score[1] / tscore) + ";q3=" + to_string(score[2] / tscore) + "]'";;
		}
		#endif
		#ifdef G_SUPPORT
		double p0 = gTest(score[0], score[1], score[2]);
		double p1 = gTest(score[1], score[2], score[0]);
		double p2 = gTest(score[2], score[0], score[1]);
		if (support == 1) res += to_string(p0);
		else {
			array<double, 3> p = {p0, p1, p2};
			res += "'[p=" + to_string(p0) + ";p2=" + to_string(p1) + ";p3=" + to_string(p2)
				+ ";s1=" + to_string(score[0]) + ";s2=" + to_string(score[1]) + ";s3=" + to_string(score[2]) + "]'";
			if (support == 3) qInfo[v] = make_tuple(score, p, get<2>(qInfo[get<0>(c)]) + "," + get<2>(qInfo[get<1>(c)]));
		}
		node->s = p0;
		#endif
		#ifdef CUSTOMIZED_LENGTH
		node->len = quad.length(score);
		if (support != 0) res += string(":") + to_string(node->len);
		#endif
		#ifdef CUSTOMIZED_ANNOTATION_LENGTH
		if (support != 0) {
			array<score_t, 5> lengths = node->annot.lengths();
			res += string(":") + to_string(lengths[0]);
		}
		#endif
		#ifdef CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH
		if (support != 0) {
			array<score_t, 5> lengths = node->annot.lengths();
			if (nodes[get<0>(c)].leafId != -1) resL += string(":") + to_string(lengths[4]);
			if (nodes[get<1>(c)].leafId != -1) resR += string(":") + to_string(lengths[3]);
			res += string(":") + to_string(lengths[0]);
			if (toplevel) {
				tempRoot = lengths[1];
				tempSibling = lengths[2];
			}
		}
		#endif
		#ifdef CASTLES
		if (support != 0) {
			CastlesNode castles(node->annot, ARG.getDoubleArg("genelength"), ARG.getIntArg("numgenetrees"));
			if (nodes[get<0>(c)].leafId != -1) resL += string(":") + to_string(castles.leftEdgeLength);
			if (nodes[get<1>(c)].leafId != -1) resR += string(":") + to_string(castles.rightEdgeLength);
			res += string(":") + to_string(castles.edgeLengthOtherwise);
			if (toplevel) {
				tempSibling = castles.siblingEdgeLength;
			}
		}
		#else
		#ifdef SUPPORT
		if (3 * score[0] > tscore) node->len = max(0.0, -log(1.5 - 1.5 * score[0] / (tscore + lambda * 2)));
		else node->len = 0;
		res += ":" + to_string(node->len);
		#endif
		#endif
		return string("(") + resL + "," + resR + ")" + res;
	}
	
	void printFreqQuadCSV(int v, int u, string top, ostream &fcsv, unordered_map<int, tuple<array<double, 3>, array<double, 3>, string> > &qInfo){
		if (nodes[v].leafId != -1) return;
		tuple<int, int, score_t> c = nodes[v].children[nodes[v].bestChild];
		printFreqQuadCSV(get<0>(c), get<1>(c), top + "," + get<2>(qInfo[u]), fcsv, qInfo);
		printFreqQuadCSV(get<1>(c), get<0>(c), top + "," + get<2>(qInfo[u]), fcsv, qInfo);
		double tscore = get<0>(get<0>(qInfo[v])) + get<1>(get<0>(qInfo[v])) + get<2>(get<0>(qInfo[v]));
		fcsv << "N" << v << "\tt1\t{" << get<2>(qInfo[get<0>(c)]) << "}|{" << get<2>(qInfo[get<1>(c)]) << "}#{" << get<2>(qInfo[u]) << "}|{" << top << "}\t"
			<< get<0>(get<1>(qInfo[v])) << "\t" << get<0>(get<0>(qInfo[v])) << "\t" << tscore << endl;
		fcsv << "N" << v << "\tt2\t{" << top << "}|{" << get<2>(qInfo[get<1>(c)]) << "}#{" << get<2>(qInfo[u]) << "}|{" << get<2>(qInfo[get<0>(c)]) << "}\t"
			<< get<1>(get<1>(qInfo[v])) << "\t" << get<1>(get<0>(qInfo[v])) << "\t" << tscore << endl;
		fcsv << "N" << v << "\tt3\t{" << get<2>(qInfo[get<0>(c)]) << "}|{" << top << "}#{" << get<2>(qInfo[u]) << "}|{" << get<2>(qInfo[get<1>(c)]) << "}\t"
			<< get<2>(get<1>(qInfo[v])) << "\t" << get<2>(get<0>(qInfo[v])) << "\t" << tscore << endl;
		
	}
	
	string printOptimalTreeWithSupport(int support, double lambda, double weight){
		#ifdef CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH
		score_t siblingL = 0, siblingR = 0, rootSum = 0, rootCount = 0;
		#else
		#ifdef CASTLES
		length_t siblingL = 0, siblingR = 0, rootSum = 0, rootCount = 0;
		#endif
		#endif
		annotTree.reset(new AnnotatedTree());
		shared_ptr<AnnotatedTree::Node> root = annotTree->addRoot();
		shared_ptr<AnnotatedTree::Node> left = annotTree->addLeft(root);
		shared_ptr<AnnotatedTree::Node> right = annotTree->addRight(root);
		right->taxon = 0;
		right->name = names[0];

		unordered_map<int, tuple<array<double, 3>, array<double, 3>, string> > qInfo;
		
		Quadrupartition quad(tripInit);
		quad.update(0, 0);
		for (size_t i = 1; i < names.size(); i++) quad.update(2, i);
		int v = hash[-taxonHash[0]];
		tuple<int, int, score_t> c = nodes[v].children[nodes[v].bestChild];
		//0|c1|c0|-
		switchSubtree(quad, get<1>(c), 2, 1);
		string resL = printOptimalSubtreeWithSupport(quad, get<0>(c), get<1>(c), support, lambda, qInfo, weight, annotTree->addLeft(left), true);
		#ifdef CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH
		if (nodes[get<0>(c)].leafId == -1) {
			siblingL = tempSibling;
			rootSum += tempRoot;
			rootCount++;
		}
		#else
		#ifdef CASTLES
		if (nodes[get<0>(c)].leafId == -1) siblingL = tempSibling;
		#endif
		#endif
		//0|c0|c1|-
		switchSubtree(quad, get<1>(c), 1, 2);
		switchSubtree(quad, get<0>(c), 2, 1);
		string resR = printOptimalSubtreeWithSupport(quad, get<1>(c), get<0>(c), support, lambda, qInfo, weight, annotTree->addRight(left), true);
		#ifdef CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH
		if (nodes[get<1>(c)].leafId == -1) {
			siblingR = tempSibling;
			rootSum += tempRoot;
			rootCount++;
		}
		if (nodes[get<0>(c)].leafId != -1) resL += string(":") + to_string(siblingR);
		if (nodes[get<1>(c)].leafId != -1) resR += string(":") + to_string(siblingL);
		string res = string(":") + to_string(rootSum/rootCount/2) + "," + names[0] + ":" + to_string(rootSum/rootCount/2);
		#else
		#ifdef CASTLES
		if (nodes[get<1>(c)].leafId == -1) siblingR = tempSibling;
		if (nodes[get<0>(c)].leafId != -1) resL += string(":") + to_string(siblingR);
		if (nodes[get<1>(c)].leafId != -1) resR += string(":") + to_string(siblingL);
		string res = string(":") + to_string(ARG.getDoubleArg("outgrouplength")/2) + "," + names[0] + ":" + to_string(ARG.getDoubleArg("outgrouplength")/2);
		#else
		string res = string(",") + names[0];
		#endif
		#endif

		if (support == 3) {
			ofstream fcsv("freqQuad.csv");
			printFreqQuadCSV(get<0>(c), get<1>(c), names[0], fcsv, qInfo);
			printFreqQuadCSV(get<1>(c), get<0>(c), names[0], fcsv, qInfo);
		}
		return string("((") + resL + "," + resR + ")" + res + ");";
	}
	*/

	void computeOptimalSubtreeWithSupport(SpeciesTree &tree, Quadrupartition& quad, int v, int u,
		unordered_map<int, tuple<array<double, 3>, array<double, 3>, string> >& qInfo, shared_ptr<SpeciesTree::Node> node, bool toplevel = false) {
		if (nodes[v].leafId != -1) {
			int i = nodes[v].leafId;
			node->taxon = i;
			node->name = names[i];
			return;
		}
		//r|u|c0c1|-
		string res;
		tuple<int, int, score_t> c = nodes[v].children[nodes[v].bestChild];
		//ru|c1|c0|-
		switchSubtree(quad, u, 1, 0);
		switchSubtree(quad, get<1>(c), 2, 1);
		computeOptimalSubtreeWithSupport(tree, quad, get<0>(c), get<1>(c), qInfo, tree.addLeft(node));
		//ru|c0|c1|-
		switchSubtree(quad, get<0>(c), 2, 1);
		switchSubtree(quad, get<1>(c), 1, 2);
		computeOptimalSubtreeWithSupport(tree, quad, get<1>(c), get<0>(c), qInfo, tree.addRight(node));
		//r|u|c1|c0
		switchSubtree(quad, u, 0, 1);
		switchSubtree(quad, get<0>(c), 1, 3);
		#ifdef CUSTOMIZED_ANNOTATION
		auto annot = quad.annotate();
		#endif
		array<double, 3> score = quad.score();
		//r|u|c0c1|-
		switchSubtree(quad, get<0>(c), 3, 2);
		double weight = tree.temp["weight"];
		score[0] /= weight; score[1] /= weight; score[2] /= weight;
		double tscore = score[0] + score[1] + score[2];

		#ifdef SUPPORT
		double support0, support1, support2;
		double lambda = tree.temp["lambda"];
		if (tscore < 1e-8) {
			support0 = 1.0 / 3; support1 = 1.0 / 3; support2 = 1.0 / 3;
		}
		else if (tscore > 100 && score[0] - max(score[1], score[2]) > 5 * sqrt(tscore)) {
			support0 = 1; support1 = 0; support2 = 0;
		}
		else if (tscore > 1000 && score[0] - max(score[1], score[2]) > 5 * sqrt((score[0] + max(score[1], score[2])) * 0.25)) {
			support0 = 1; support1 = 0; support2 = 0;
		}
		else {
			double i0 = 1.0 - incbeta(score[0] + 1.0, tscore + lambda * 2 - score[0], 1.0 / 3.0);
			double i1 = 1.0 - incbeta(score[1] + 1.0, tscore + lambda * 2 - score[1], 1.0 / 3.0);
			double i2 = 1.0 - incbeta(score[2] + 1.0, tscore + lambda * 2 - score[2], 1.0 / 3.0);
			double lb0 = lgamma(score[0] + 1.0) + lgamma(tscore - score[0] + lambda * 2);
			double lb1 = lgamma(score[1] + 1.0) + lgamma(tscore - score[1] + lambda * 2);
			double lb2 = lgamma(score[2] + 1.0) + lgamma(tscore - score[2] + lambda * 2);
			support0 = i0 / (i0 + i1 * exp(log(2.0) * (score[1] - score[0]) + lb1 - lb0) + i2 * exp(log(2.0) * (score[2] - score[0]) + lb2 - lb0));
			support1 = i1 / (i1 + i0 * exp(log(2.0) * (score[0] - score[1]) + lb0 - lb1) + i2 * exp(log(2.0) * (score[2] - score[1]) + lb2 - lb1));
			support2 = i2 / (i2 + i1 * exp(log(2.0) * (score[1] - score[2]) + lb1 - lb2) + i0 * exp(log(2.0) * (score[0] - score[2]) + lb0 - lb2));
		}
		node->attributes["support"] = support0;
		node->attributes["pp1"] = support0;
		node->attributes["pp2"] = support1;
		node->attributes["pp3"] = support2;
		node->attributes["f1"] = score[0];
		node->attributes["f2"] = score[1];
		node->attributes["f3"] = score[2];
		node->attributes["q1"] = score[0] / tscore;
		node->attributes["q2"] = score[1] / tscore;
		node->attributes["q3"] = score[2] / tscore;
		//if (support == 3) qInfo[v] = make_tuple(score, p, get<2>(qInfo[get<0>(c)]) + "," + get<2>(qInfo[get<1>(c)]));
		#endif
		#ifdef LOCAL_BOOTSTRAP
		array<int, 3> bs = annot.bootstrap();
		double nbs = annot.bs.size() / 100.0;
		node->attributes["support"] = bs[0] / nbs;
		node->attributes["bs1"] = bs[0];
		node->attributes["bs2"] = bs[1];
		node->attributes["bs3"] = bs[2];
		node->attributes["s1"] = score[0];
		node->attributes["s2"] = score[1];
		node->attributes["s3"] = score[2];
		node->attributes["q1"] = score[0] / tscore;
		node->attributes["q2"] = score[1] / tscore;
		node->attributes["q3"] = score[2] / tscore;
		#endif
		#ifdef G_SUPPORT
		double p0 = gTest(score[0], score[1], score[2]);
		double p1 = gTest(score[1], score[2], score[0]);
		double p2 = gTest(score[2], score[0], score[1]);
		if (support == 1) res += to_string(p0);
		else {
			array<double, 3> p = { p0, p1, p2 };
			res += "'[p=" + to_string(p0) + ";p2=" + to_string(p1) + ";p3=" + to_string(p2)
				+ ";s1=" + to_string(score[0]) + ";s2=" + to_string(score[1]) + ";s3=" + to_string(score[2]) + "]'";
			if (support == 3) qInfo[v] = make_tuple(score, p, get<2>(qInfo[get<0>(c)]) + "," + get<2>(qInfo[get<1>(c)]));
		}
		node->s = p0;
		#endif
		#ifdef CUSTOMIZED_LENGTH
		node->attributes["length"] = quad.length(score);
		#endif
		#ifdef CUSTOMIZED_ANNOTATION_LENGTH
		array<score_t, 5> lengths = annot.lengths();
		node->attributes["length"] = lengths[0];
		#endif
		#ifdef CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH
		array<score_t, 5> lengths = annot.lengths();
		node->attributes["length"] = lengths[0];
		if (nodes[get<0>(c)].leafId != -1) node->lc->attributes["length"] = lengths[4];
		if (nodes[get<1>(c)].leafId != -1) node->rc->attributes["length"] = lengths[3];
		if (toplevel) {
			tree.temp["rootlengthsum"] += lengths[1];
			tree.temp["rootlengthcnt"]++;
			tree.temp["rootsiblinglength"] = lengths[2];
		}
		#endif
		#ifdef CASTLES
		CastlesNode castles(annot, ARG.getDoubleArg("genelength"), ARG.getIntArg("numgenetrees"));
		node->attributes["length"] = castles.edgeLengthOtherwise;
		if (nodes[get<0>(c)].leafId != -1) node->lc->attributes["length"] = castles.leftEdgeLength;
		if (nodes[get<1>(c)].leafId != -1) node->rc->attributes["length"] = castles.rightEdgeLength;
		if (toplevel) {
			tree.temp["rootsiblinglength"] = castles.siblingEdgeLength;
		}
		#else
		#ifdef SUPPORT
		if (3 * score[0] > tscore) node->attributes["length"] = max(0.0, -log(1.5 - 1.5 * score[0] / (tscore + lambda * 2)));
		else node->attributes["length"] = 0;
		#endif
		#endif
	}

	SpeciesTree optimalTreeWithSupport(double lambda, double weight) {
		SpeciesTree tree;
		tree.temp["lambda"] = lambda;
		tree.temp["weight"] = weight;

		#if defined(CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH) || defined(CASTLES)
		double siblingL = 0, siblingR = 0, rootSum = 0, rootCount = 0;
		#endif
		shared_ptr<SpeciesTree::Node> root = tree.addRoot();
		shared_ptr<SpeciesTree::Node> left = tree.addLeft(root);
		shared_ptr<SpeciesTree::Node> right = tree.addRight(root);
		right->taxon = 0;
		right->name = names[0];

		unordered_map<int, tuple<array<double, 3>, array<double, 3>, string> > qInfo;

		Quadrupartition quad(tripInit);
		quad.update(0, 0);
		for (size_t i = 1; i < names.size(); i++) quad.update(2, i);
		int v = hash[-taxonHash[0]];
		tuple<int, int, score_t> c = nodes[v].children[nodes[v].bestChild];
		//0|c1|c0|-
		switchSubtree(quad, get<1>(c), 2, 1);
		computeOptimalSubtreeWithSupport(tree, quad, get<0>(c), get<1>(c), qInfo, tree.addLeft(left), true);
		//0|c0|c1|-
		switchSubtree(quad, get<1>(c), 1, 2);
		switchSubtree(quad, get<0>(c), 2, 1);
		computeOptimalSubtreeWithSupport(tree, quad, get<1>(c), get<0>(c), qInfo, tree.addRight(left), true);
		#if defined(CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH) || defined(CASTLES)
		if (nodes[get<0>(c)].leafId != -1) left->lc->attributes["length"] = tree.temp["rootsiblinglength"];
		if (nodes[get<1>(c)].leafId != -1) left->rc->attributes["length"] = tree.temp["rootsiblinglength"];
		#ifdef CASTLES
		left->attributes["length"] = ARG.getDoubleArg("outgrouplength");
		#else
		left->attributes["length"] = tree.temp["rootlengthsum"] / tree.temp["rootlengthcnt"];
		#endif
		#endif

		/*
		if (support == 3) {
			ofstream fcsv("freqQuad.csv");
			printFreqQuadCSV(get<0>(c), get<1>(c), names[0], fcsv, qInfo);
			printFreqQuadCSV(get<1>(c), get<0>(c), names[0], fcsv, qInfo);
		}
		*/

		return tree;
	}

	string subtreeLeaves(shared_ptr<SpeciesTree::Node> node){
		if (node->isLeaf()) return node->name;
		else return node->name = subtreeLeaves(node->lc) + "," + subtreeLeaves(node->rc);
	}
	
	void printFreqQuadCSVRecursion(shared_ptr<SpeciesTree::Node> v, shared_ptr<SpeciesTree::Node> u, string top, ostream &fcsv, int &id){
		if (v->isLeaf()) return;
		printFreqQuadCSVRecursion(v->lc, v->rc, top + "," + u->name, fcsv, id);
		printFreqQuadCSVRecursion(v->rc, v->lc, top + "," + u->name, fcsv, id);
		id++;
		double tscore = v->attributes["f1"] + v->attributes["f2"] + v->attributes["f3"]; 
		fcsv << "N" << id << "\tt1\t{" << v->lc->name << "}|{" << v->rc->name << "}#{" << u->name << "}|{" << top << "}\t"
			<< v->attributes["pp1"] << "\t" << v->attributes["f1"] << "\t" << tscore << endl;
		fcsv << "N" << id << "\tt2\t{" << top << "}|{" << v->rc->name << "}#{" << u->name << "}|{" << v->lc->name << "}\t"
			<< v->attributes["pp2"] << "\t" << v->attributes["f2"] << "\t" << tscore << endl;
		fcsv << "N" << id << "\tt3\t{" << v->lc->name << "}|{" << top << "}#{" << u->name << "}|{" << v->rc->name << "}\t"
			<< v->attributes["pp3"] << "\t" << v->attributes["f3"] << "\t" << tscore << endl;
	}
	
	void printFreqQuadCSV(SpeciesTree& tree){
		subtreeLeaves(tree.root());
		ofstream fcsv("freqQuad.csv");
		int id = 0;
		printFreqQuadCSVRecursion(tree.root()->lc->lc, tree.root()->lc->rc, tree.root()->rc->name, fcsv, id);
		printFreqQuadCSVRecursion(tree.root()->lc->rc, tree.root()->lc->lc, tree.root()->rc->name, fcsv, id);
	}

#endif
};

struct MetaAlgorithm{
	static bool STATIC_INITIALIZED;

	vector<string> files, names;
	int nThreads = 1, nRounds = 4, nSample = 4, support = 1;
	double p = 0.25, lambda = 0.5;
	string outputFile, guideFile, constraintFile, constraintTree, guideTree;
	ofstream fileOut;
	unordered_map<string, int> name2id;
	unordered_map<string, string> leafname_mapping;
	
	TripartitionInitializer tripInit;

	shared_ptr<AnnotatedTree> annotTree;
	
	MetaAlgorithm(){
		initialize();
	}
	
	MetaAlgorithm(int argc, char** argv) {
		staticInitialize(argc, argv);
		initialize();
	}
	
	static void staticInitialize(int argc, char** argv) {
		if (STATIC_INITIALIZED) return;
		STATIC_INITIALIZED = true;

		string version = string(ALG_VERSION) + "." + OBJECTIVE_VERSION + "." + DRIVER_VERSION;
		MDGenerator::version = version;

		cerr << ARG.getFullName() << endl;
		cerr << "Version: " << version << endl;
		ARG.addStringArg('c', "constraint", "", "Newick file containing a binary species tree to place missing species on");
		ARG.addStringArg('g', "guide", "", "Newick file containing binary trees as guide trees");
		ARG.addStringArg('o', "output", "<standard output>", "File name for the output species tree", true);
		ARG.addIntArg('r', "round", 4, "Number of initial rounds of placements");
		ARG.addIntArg('s', "subsample", 4, "Number of rounds of subsampling per exploration step");
		ARG.addDoubleArg(0, "proportion", 0.25, "Proportion of taxa in the subsample in naive algorithm");
		ARG.addIntArg('t', "thread", 1, "Number of threads", true);
		ARG.addIntArg(0, "seed", 233, "Seed for pseudorandomness");
		ARG.addIntArg('v', "verbose", 2, "Level of logging (1: minimum, 2: normal)");
		ARG.addFlag('C', "scoring", "Scoring the full species tree file after `-c` without exploring other topologies (`-r 1 -s 0`)", [&]() {
			ARG.getIntArg("round") = 1; ARG.getIntArg("subsample") = 0;
			}, true);
		ARG.addFlag('R', "moreround", "More rounds of placements and subsampling (`-r 16 -s 16`)", [&]() {
			ARG.getIntArg("round") = 16; ARG.getIntArg("subsample") = 16;
			}, true);
	#ifdef SUPPORT
		ARG.addDoubleArg('l', "lambda", 0.5, "Rate lambda of Yule process under which the species tree is modeled");
		ARG.addDoubleArg('w', "downweightrepeat", 1, "The number of trees sampled for each locus");
		ARG.addIntArg('u', "support", 1, "Output support option (0: no output support value, 1: branch local posterior probability, 2: detailed, 3: freqQuad.csv)");
	#endif
	#ifdef G_SUPPORT
		ARG.addIntArg('u', "support", 1, "Output support option (0: no output support value, 1: Grubbs's test p-value, 2: detailed, 3: freqQuad.csv)");
	#endif
	#ifdef LOCAL_BOOTSTRAP
		ARG.addIntArg('u', "support", 1, "Output support option (0: no output support value, 1: local bootstrap, 2: detailed)");
	#endif
	#ifdef ROOTING
		ARG.addStringArg(0, "root", "", "Root at the given species");
	#endif
	#ifdef NAME_MAPPING
		ARG.addStringArg('a', "mapping", "", "A list of gene name to taxon name maps, each line contains one gene name followed by one taxon name separated by a space or tab");
	#endif
	
		ARG.parse(argc, argv);
	}

	void initialize(int argc, char** argv, string dummy1 = "", string dummy2 = "") {
		staticInitialize(argc, argv);
		initialize();
	}

	void initialize(){
		constraintFile = ARG.getStringArg("constraint");
		guideFile = ARG.getStringArg("guide");
		outputFile = ARG.getStringArg("output");
		nRounds = ARG.getIntArg("round");
		nSample = ARG.getIntArg("subsample");
		p = ARG.getDoubleArg("proportion");
		nThreads = ARG.getIntArg("thread");
		lambda = ARG.getDoubleArg("lambda");
		support = ARG.getIntArg("support");
		srand(ARG.getIntArg("seed"));
		RND_GENERATOR.seed(ARG.getIntArg("seed"));

		int loglevel = ARG.getIntArg("verbose");
		LOG.enabled = (loglevel >= 2);
	
	#ifdef NAME_MAPPING
		string mappingFile;
		mappingFile = ARG.getStringArg("mapping");
		if (mappingFile != ""){
			ifstream fmap(mappingFile);
			string gname, sname;
			while (fmap >> gname){
				fmap >> sname;
				leafname_mapping[gname] = sname;
			}
		}
	#endif

	}

	string mappedname(string name){
		if (leafname_mapping.count(name)) return leafname_mapping.at(name);
		return name;
	}
	
	pair<score_t, string> run(){
		ThreadPool TP;
		TP.initialize(nThreads);
		int ROUND_NN = -1;

		ostream &fout = (outputFile == "" || outputFile == "<standard output>") ? cout : fileOut;
		if (outputFile != "" && outputFile != "<standard output>") fileOut.open(outputFile);
		// fout << setprecision(3);

		LOG << "#Species: " << names.size() << endl;
		LOG << "#Rounds: " << nRounds << endl;
		LOG << "#Samples: " << nSample << endl;
		LOG << "#Threads: " << nThreads << endl;
		
		if (constraintFile != ""){
			ifstream fin(constraintFile);
			getline(fin, constraintTree);
		}
		else{
			ROUND_NN = 20 + 2 * sqrt(names.size()) * log2(names.size());
		}
		
		ConstrainedOptimizationAlgorithm alg(names.size(), tripInit, names, TP, ROUND_NN);
		
		if (guideFile != ""){
			ifstream fin(guideFile);
			string tree;
			while (getline(fin, tree)){
				if (tree.size()) alg.addGuideTree(tree, name2id, ROUND_NN);
			}
		}
		if (guideTree != ""){
			alg.addGuideTree(guideTree, name2id, ROUND_NN);
		}
		
		auto res = (constraintTree == "") ? alg.run(nRounds) : alg.constrainedRun(nRounds, constraintTree, name2id);
		LOG << "Initial score: " << (double) res.first << endl;
		LOG << "Initial tree: " << res.second << endl;
		if (constraintTree == "") {
			LOG << "*** Subsample Process ***" << endl;
			score_t prevS;
			int roundNum = 0;
			do {
				prevS = res.first;
				res = alg.run(nSample, p);
				LOG << "Current score: " << (double) res.first << endl;
				LOG << "Current tree: " << res.second << endl;
				roundNum++;
			}
			while (prevS + ERROR_TOLERANCE < res.first && roundNum < 20);
			if (prevS + ERROR_TOLERANCE < res.first) LOG << "Search stopped due to excessive rounds. Tree may not be optimal!" << endl;
		}
		
		string output = res.second;
		LOG << "Final Tree: " << output << endl;
		#ifdef SOME_FORM_OF_SUPPORT
		double w = 1;
		#ifdef SUPPORT
		w = ARG.getDoubleArg("downweightrepeat");
		#endif
		if (support) {
			SpeciesTree tree = alg.optimalTreeWithSupport(lambda, w);
			if (support == 1) output = tree.simpleTree("length", "support");
			else output = tree.annotatedTree("length");
			if (support == 3) alg.printFreqQuadCSV(tree);
		}
		#endif
		fout << output << endl;

		annotTree = alg.annotTree;

		return res;
	}
};

bool MetaAlgorithm::STATIC_INITIALIZED = false;
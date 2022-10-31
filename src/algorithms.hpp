#define ALG_VERSION "v1.10"

/* CHANGE LOG
 * 1.10: fixing constraint/guide trees for rooting on nodes
 * 1.9: deleting obsolete code
 * 1.8: adding seed to pseudorandomness
 * 1.7: adding Grubbs's test based support
 * 1.6: adding -w option for gene tree replications
 * 1.5: improving parallelization
 */

typedef unsigned __int128 hash_t;

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
#include<mutex>
#include<cmath>

#ifdef SUPPORT
#include "incbeta.c"
#endif

#ifndef ERROR_TOLERANCE
#define ERROR_TOLERANCE 0
#endif
using namespace std;

ThreadPool TP;

struct Hasher{
	size_t operator()(const hash_t h) const{
		return (size_t) h;
	}
};

int ROUND_NN = -1;

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
	int rootNodeId = -1, rootLeafId = -1, orderId = 0;
	int rNN;
	
	Tripartition trip;
	
	PlacementAlgorithm(const vector<hash_t> &taxonHash, const TripartitionInitializer &tripInit, int rNN = ROUND_NN): taxonHash(taxonHash), trip(tripInit), rNN(rNN){}
	
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
			//trip.addTotal(order[orderId]);
			orderId++;
		}
		if (rootNodeId == -1){
			rootNodeId = nodes.size();
			nodes.emplace_back(order[orderId]);
			//trip.addTotal(order[orderId]);
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
		
		//E|AB|CD
		switchSubtree(light(v), 1, 2);
		score_t abe_c_d = 0, cde_a_b = 0;
		if (leafId(light(v)) == -1){
			int c = heavy(light(v)), d = light(light(v));
			//ABE|C|D
			switchSubtreeSet(heavy(v), 1, 0);
			switchSubtreeSet(c, 2, 1);
			tripScoreSet();
			//ABE|CD|-
			switchSubtreeSet(d, 2, 1);

			//ABE|C|D
			switchSubtreeGet(heavy(v), 1, 0);
			switchSubtreeGet(c, 2, 1);
			abe_c_d = tripScoreGet();
			//ABE|CD|-
			switchSubtreeGet(d, 2, 1);
			if (recurse) abe_c_d = nnLocal(light(v), false, abe_c_d);
			//E|AB|CD
			switchSubtreeSet(heavy(v), 0, 1);
			switchSubtreeSet(light(v), 1, 2);

			//E|AB|CD
			switchSubtreeGet(heavy(v), 0, 1);
			switchSubtreeGet(light(v), 1, 2);
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
			//ABE|CD|-
			score_t abe_c_d = nnMove(light(v));
			score_t ae_b_cd = 0;
			if (leafId(heavy(v)) == -1){
				int b = light(heavy(v));
				//AE|B|CD
				switchSubtreeSet(b, 0, 1);
				switchSubtreeSet(light(v), 1, 2);
				tripScoreSet();
				//ABCDE|-|-
				switchSubtreeSet(b, 1, 0);
				switchSubtreeSet(light(v), 2, 0);

				//AE|B|CD
				switchSubtreeGet(b, 0, 1);
				switchSubtreeGet(light(v), 1, 2);
				ae_b_cd = tripScoreGet();
				//ABCDE|-|-
				switchSubtreeGet(b, 1, 0);
				switchSubtreeGet(light(v), 2, 0);
			}
			//ABCDE|-|-
			else switchSubtree(light(v), 1, 0);
			//CDE|AB|-
			score_t cde_a_b = nnMove(heavy(v));
			//E|AB|CD
			switchSubtree(light(v), 0, 2);
			score_t e_ab_cd = tripScore();
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
				int b = light(heavy(v));
				//E|A|BCD
				switchSubtreeSet(b, 1, 2);
				tripScoreSet();
				//BE|A|CD
				switchSubtreeSet(b, 2, 0);
				tripScoreSet();
				//E|ACD|B
				switchSubtreeSet(b, 0, 2);
				switchSubtreeSet(light(v), 2, 1);
				tripScoreSet();
				//E|AB|CD
				switchSubtreeSet(b, 2, 1);
				switchSubtreeSet(light(v), 1, 2);

				//E|A|BCD
				switchSubtreeGet(b, 1, 2);
				score_t e_a_bcd = tripScoreGet();
				//BE|A|CD
				switchSubtreeGet(b, 2, 0);
				score_t be_a_cd = tripScoreGet();
				//E|ACD|B
				switchSubtreeGet(b, 0, 2);
				switchSubtreeGet(light(v), 2, 1);
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
				//E|AB|CD
				switchSubtreeGet(b, 2, 1);
				switchSubtreeGet(light(v), 1, 2);
			}
			//E|ABCD|-
			switchSubtree(light(v), 2, 1);
			return nnPerform(v, best_case, top_score);
		}
	}
	
	void nnMove(){
		if (ROUND_NN >= 0) {
			tripUpdateSet(0, rootLeafId);
			switchSubtreeSet(rootNodeId, -1, 0);
			tripUpdateGet(0, rootLeafId);
			switchSubtreeGet(rootNodeId, -1, 0);
			nnMove(rootNodeId);
			cerr << "#NNI moves:" << ROUND_NN - rNN << "/" << ROUND_NN << endl;
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
			if (orderId & 15 == 0) cerr << "Placing " << orderId << "/" << order.size() << endl;
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
		Node(int leafId, hash_t hash): leafId(leafId), hash(hash){}
		
		hash_t hash, bestChild = 0;
		score_t bestScore = 0;
		int roundId = 0, leafId;
		unordered_map<hash_t, tuple<int, int, score_t>, Hasher> children;
	};
	
	vector<Node> nodes;
	unordered_map<hash_t, int, Hasher> hash;
	vector<hash_t> taxonHash;
	const vector<string> names;
	mt19937_64 generator;
	uniform_int_distribution<hash_t> randomHash;
	uniform_real_distribution<double> randP;
	const TripartitionInitializer &tripInit;
	const int ntaxa;
	int roundId = 0;
	mutex mtx;
	
	ConstrainedOptimizationAlgorithm(const int ntaxa, const TripartitionInitializer &tripInit, const vector<string> &names, const int seed = rand()): 
			ntaxa(ntaxa), tripInit(tripInit), names(names), generator(seed){
		taxonHash.push_back(0);
		for (int i = 1; i < ntaxa; i++){
			hash_t r = randomHash(generator);
			taxonHash.push_back(r);
			hash[r] = nodes.size();
			nodes.emplace_back(i, r);
			taxonHash[0] -= r;
		}
	}
	
	ConstrainedOptimizationAlgorithm(const ConstrainedOptimizationAlgorithm &alg): 
			ntaxa(alg.ntaxa), tripInit(alg.tripInit), names(alg.names), nodes(alg.nodes), hash(alg.hash), taxonHash(alg.taxonHash), generator(rand()){}
	
	int subsampleSubtree(int v, PlacementAlgorithm &pAlg, const unordered_set<int> &selected){
		if (nodes[v].leafId != -1){
			if (selected.count(nodes[v].leafId)){
				//pAlg.trip.addTotal(nodes[v].leafId);
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
				//pAlg.trip.addTotal(0);
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
		shuffle(pAlg.order.begin(), pAlg.order.end(), generator);
	}
	
	int subsampleSubtree(int v, PlacementAlgorithm &pAlg, double subsampleRate){
		if (nodes[v].leafId != -1){
			if (randP(generator) < subsampleRate){
				//pAlg.trip.addTotal(nodes[v].leafId);
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
			if (randP(generator) < subsampleRate){
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
		shuffle(pAlg.order.begin(), pAlg.order.end(), generator);
	}
	
	PlacementAlgorithm createPlacementAlgorithm(double subsampleRate){
		PlacementAlgorithm pAlg(taxonHash, tripInit);
		createPlacementAlgorithm(pAlg, subsampleRate);
		return pAlg;
	}
	
	int guideSubtree(PlacementAlgorithm &pAlg, const string &tree, const unordered_map<string, int> &name2id, int &i, unordered_set<int> &added){
		int ret;
		if (tree[i] != '('){
			string s;
			for (; tree[i] != ':' && tree[i] != ')' && tree[i] != ',' && tree[i] != ';'; i++){
				if (tree[i] != '\"' && tree[i] != '\'') s += tree[i];
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
				cerr << "Warning: Polytomy detected! Resolving arbitrarily!\n";
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
	
	void addGuideTree(const string tree, const unordered_map<string, int> &name2id, int rnn = ROUND_NN){
		unordered_set<int> added;
		PlacementAlgorithm pAlg(taxonHash, tripInit, rnn);
		int i = 0;
		guideSubtree(pAlg, tree, name2id, i, added);
		for (int i = 0; i < ntaxa; i++){
			if (added.count(i) == 0) pAlg.order.push_back(i);
		}
		{ const lock_guard<mutex> lock(mtx); shuffle(pAlg.order.begin(), pAlg.order.end(), generator); }
		pAlg.run();
		cerr << pAlg.printTree(names) << endl;
		{ const lock_guard<mutex> lock(mtx); addTripartitions(pAlg.tripHash); }
	}
	
	int hash2node(hash_t h){
		if (hash.count(h) == 0){
			hash[h] = nodes.size();
			nodes.emplace_back(-1, h);
		}
		return hash[h];
	}
	
	void addTripartitions(const vector<tuple<hash_t, hash_t, score_t> > tripHash){
		for (const tuple<hash_t, hash_t, score_t> t: tripHash){
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
	
	void batchWork(vector<PlacementAlgorithm> &jobs, int start, int end){
		for (int i = start; i < end; i++){
			jobs[i].run();
			cerr << jobs[i].printTree(names) << endl;
		}
	}
	
	pair<score_t, string> run(int nJobs = 1, int nThrds = 1, double subsampleRate = 0, bool allowTwoStepRun = true){
		if (nJobs == 0) return {computeOptimalTree(), printOptimalTree()};
		if (allowTwoStepRun && ntaxa >= 100) return twoStepRun(nJobs, nThrds);
		vector<PlacementAlgorithm> jobs;
		vector<thread> thrds;
		for (int i = 0; i < nJobs; i++){
			jobs.emplace_back(move(createPlacementAlgorithm(subsampleRate)));
		}
		for (int i = 1; i < nThrds; i++){
			thrds.emplace_back(&ConstrainedOptimizationAlgorithm::batchWork, this, ref(jobs), i * nJobs / nThrds, (i + 1) * nJobs / nThrds);
		}
		batchWork(jobs, 0, nJobs / nThrds);
		for (int i = 1; i < nThrds; i++){
			thrds[i - 1].join();
		}
		for (int i = 0; i < nJobs; i++){
			addTripartitions(jobs[i].tripHash);
		}
		
		return {computeOptimalTree(), printOptimalTree()};
	}
	
	void batchConstrainedWork(int start, int end, const string tree, const unordered_map<string, int> &name2id){
		for (int i = start; i < end; i++){
			addGuideTree(tree, name2id, -1);
		}
	}
	
	pair<score_t, string> constrainedRun(int nJobs, int nThrds, const string tree, const unordered_map<string, int> &name2id){
		vector<thread> thrds;
		for (int i = 1; i < nThrds; i++){
			thrds.emplace_back(&ConstrainedOptimizationAlgorithm::batchConstrainedWork, this, i * nJobs / nThrds, (i + 1) * nJobs / nThrds, tree, ref(name2id));
		}
		batchConstrainedWork(0, nJobs / nThrds, tree, name2id);
		for (int i = 1; i < nThrds; i++){
			thrds[i - 1].join();
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
			//pAlg.trip.update(-1, i);
			pAlg.tripUpdateSet(-1, i);
			pAlg.tripUpdateGet(-1, i);
			return false;
		}
	}
	
	void twoStepWorkflow(PlacementAlgorithm &pAlg){
		int N, rId;
		{ const lock_guard<mutex> lock(mtx); N = ntaxa; rId = roundId; }
		int n = sqrt(N) * log2(names.size()) / 4;
		vector<int> order;
		if (rId == 0){
			mtx.lock();
			ConstrainedOptimizationAlgorithm alg(*this);
			for (int i = 1; i < N; i++) order.push_back(i);
			shuffle(order.begin(), order.end(), generator);
			order.push_back(order[0]);
			order[0] = 0;
			hash_t hashsum = 0;
			for (int i = 1; i < n; i++) hashsum += taxonHash[order[i]];
			mtx.unlock();
			
			pAlg.taxonHash[0] = -hashsum;
			alg.taxonHash[0] = -hashsum;
			for (int r = 0; r < n; r++){
				cerr << "Guide Tree " << r << "/" << n << endl;
				for (int i = 0; i < n; i++) pAlg.order.push_back(order[i]);
				{ const lock_guard<mutex> lock(mtx); shuffle(pAlg.order.begin(), pAlg.order.end(), generator); }
				pAlg.run();
				alg.addTripartitions(pAlg.tripHash);
				//for (int i = 0; i < n; i++) pAlg.trip.update(-1, order[i]);
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
			//cerr << alg.printOptimalTree() << endl;
			alg.createPlacementAlgorithm(pAlg, 1);
			{ const lock_guard<mutex> lock(mtx); pAlg.taxonHash[0] = taxonHash[0]; }
		}
		else{
			for (int i = 0; i < N; i++) order.push_back(i);
			{ const lock_guard<mutex> lock(mtx); shuffle(order.begin(), order.end(), generator); }
			unordered_set<int> selected;
			for (int i = 0; i < n; i++) selected.insert(order[i]);
			{ const lock_guard<mutex> lock(mtx); createPlacementAlgorithm(pAlg, selected); }
			pAlg.nnMove();
		}
		{
			vector<int> abnormalOrder;
			vector<PlacementAlgorithm::Node> pNodes(pAlg.nodes);
			int rootNodeId = pAlg.rootNodeId;//, rootLeafId = pAlg.rootLeafId;
			{
				vector<PlacementAlgorithm::Node> originalNodes(pAlg.nodes);
				int originalRootNodeId = pAlg.rootNodeId;//, originalRootLeafId = pAlg.rootLeafId;
				
				vector<vector<int> > nodeBranch(pNodes.size());
				
				for (int i = n; i < N; i++) {
					if ((i - n) & 127 == 0) cerr << "Binning " << i - n << "/" << N - n << endl; 
					nodeBranch[pAlg.locateBranch(order[i])].push_back(order[i]);
				}
				int cnt = 0;
				for (int b = 0; b < nodeBranch.size(); b++){
					pAlg.nodes = originalNodes;
					pAlg.rootNodeId = originalRootNodeId;
					//pAlg.rootLeafId = originalRootLeafId;
					vector<int> added;
					unordered_map<int, int> nodeRemap;
					nodeRemap[b] = b;
					for (int i: nodeBranch[b]) {
						if (cnt & 127 == 0) cerr << "Placement " << cnt++ << "/" << N - n << endl;
						if (placeNode(pAlg, pNodes, rootNodeId, i, nodeRemap)) added.push_back(i);
						else abnormalOrder.push_back(i);
					}
					//for (int i: added) pAlg.trip.update(-1, i);
					for (int i: added) pAlg.tripUpdateSet(-1, i);
					for (int i: added) pAlg.tripUpdateGet(-1, i);
				}
			}
			pAlg.order = abnormalOrder;
			pAlg.nodes = pNodes;
			pAlg.rootNodeId = rootNodeId;
			//pAlg.rootLeafId = rootLeafId;
			pAlg.orderId = 0;
		}
		cerr << "Remaining: " << pAlg.order.size() << endl;
		pAlg.nnMove();
		pAlg.run();
	}
	
	void twoStepBatchWork(int nJobs){
		for (int i = 0; i < nJobs; i++){
			PlacementAlgorithm job(taxonHash, tripInit);
			twoStepWorkflow(job);
			cerr << job.printTree(names) << endl;
			{ const lock_guard<mutex> lock(mtx); addTripartitions(job.tripHash); }
		}
	}
	
	pair<score_t, string> twoStepRun(int nJobs, int nThrds){
		cerr << "Use two-step algorithm!\n";
		vector<thread> thrds;
		
		for (int i = 1; i < nThrds; i++){
			thrds.emplace_back(&ConstrainedOptimizationAlgorithm::twoStepBatchWork, this, (i + 1) * nJobs / nThrds - i * nJobs / nThrds);
		}
		twoStepBatchWork(nJobs / nThrds);
		for (int i = 1; i < nThrds; i++){
			thrds[i - 1].join();
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

#if defined(SUPPORT) || defined(G_SUPPORT)
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

	string printOptimalSubtreeWithSupport(Quadrupartition &quad, int v, int u, int support, double lambda,
			unordered_map<int, tuple<array<double, 3>, array<double, 3>, string> > &qInfo, double weight){
		if (nodes[v].leafId != -1) {
			array<double, 3> t;
			if (support == 3) qInfo[v] = make_tuple(t, t, names[nodes[v].leafId]);
			return names[nodes[v].leafId];
		}
		//r|u|c0c1|-
		string res = "(";
		tuple<int, int, score_t> c = nodes[v].children[nodes[v].bestChild];
		//ru|c1|c0|-
		switchSubtree(quad, u, 1, 0);
		switchSubtree(quad, get<1>(c), 2, 1);
		res += printOptimalSubtreeWithSupport(quad, get<0>(c), get<1>(c), support, lambda, qInfo, weight) + ",";
		//ru|c0|c1|-
		switchSubtree(quad, get<0>(c), 2, 1);
		switchSubtree(quad, get<1>(c), 1, 2);
		res += printOptimalSubtreeWithSupport(quad, get<1>(c), get<0>(c), support, lambda, qInfo, weight) + ")";
		//r|u|c1|c0
		switchSubtree(quad, u, 0, 1);
		switchSubtree(quad, get<0>(c), 1, 3);
		array<double, 3> score = quad.score();
		//r|u|c0c1|-
		switchSubtree(quad, get<0>(c), 3, 2);
		score[0] /= weight; score[1] /= weight; score[2] /= weight;
		double tscore = score[0] + score[1] + score[2];
		double support0, support1, support2;
		
		#ifdef SUPPORT
		if (tscore < 1e-8){
			support0 = 1.0 / 3; support1 = 1.0 / 3; support2 = 1.0 / 3;
		}
		else if (tscore > 100 && score[0] - max(score[1], score[2]) > 5 * sqrt(tscore)){
			support0 = 1; support1 = 0; support2 = 0;
		}
		else {
			double i0 = 1.0 - incbeta(score[0] + 1.0, tscore + lambda * 2 - score[0], 1.0 / 3.0);
			double i1 = 1.0 - incbeta(score[1] + 1.0, tscore + lambda * 2 - score[1], 1.0 / 3.0);
			double i2 = 1.0 - incbeta(score[2] + 1.0, tscore + lambda * 2 - score[2], 1.0 / 3.0);
			double lb0 = lgamma(score[0] + 1.0) + lgamma(tscore - score[0] + lambda * 2) - lgamma(tscore + 1.0 + lambda * 2);
			double lb1 = lgamma(score[1] + 1.0) + lgamma(tscore - score[1] + lambda * 2) - lgamma(tscore + 1.0 + lambda * 2);
			double lb2 = lgamma(score[2] + 1.0) + lgamma(tscore - score[2] + lambda * 2) - lgamma(tscore + 1.0 + lambda * 2);
			support0 = i0 / (i0 + i1 * exp(log(2.0) * (score[1] - score[0]) + lb1 - lb0) + i2 * exp(log(2.0) * (score[2] - score[0]) + lb2 - lb0));
			support1 = i1 / (i1 + i0 * exp(log(2.0) * (score[0] - score[1]) + lb0 - lb1) + i2 * exp(log(2.0) * (score[2] - score[1]) + lb2 - lb1));
			support2 = i2 / (i2 + i1 * exp(log(2.0) * (score[1] - score[2]) + lb1 - lb2) + i0 * exp(log(2.0) * (score[0] - score[2]) + lb0 - lb2));
		}
		if (support == 1) res += to_string(support0);
		else {
			array<double, 3> p = {support0, support1, support2};
			res += "'[pp1=" + to_string(p[0]) + ";pp2=" + to_string(p[1]) + ";pp3=" + to_string(p[2]) + ";f1=" + to_string(score[0]) + ";f2=" + to_string(score[1]) + ";f3=" + to_string(score[2]) + "]'";
			if (support == 3) qInfo[v] = make_tuple(score, p, get<2>(qInfo[get<0>(c)]) + "," + get<2>(qInfo[get<1>(c)]));
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
		#endif
		#ifdef CUSTOMIZED_LENGTH
		if (support != 0){
			res += ":";
			res += quad.length(score);
		}
		#endif
		#if defined(SUPPORT)
		if (3 * score[0] > tscore) res += ":" + to_string(max(0.0, -log(1.5 - 1.5 * score[0] / (tscore + lambda * 2))));
		else res += ":" + to_string(0.0);
		#endif
		return res;
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
		unordered_map<int, tuple<array<double, 3>, array<double, 3>, string> > qInfo;
		string res;
		Quadrupartition quad(tripInit);
		quad.update(0, 0);
		for (int i = 1; i < names.size(); i++) quad.update(2, i);
		int v = hash[-taxonHash[0]];
		tuple<int, int, score_t> c = nodes[v].children[nodes[v].bestChild];
		//0|c1|c0|-
		switchSubtree(quad, get<1>(c), 2, 1);
		res = "((" + printOptimalSubtreeWithSupport(quad, get<0>(c), get<1>(c), support, lambda, qInfo, weight);
		//0|c0|c1|-
		switchSubtree(quad, get<1>(c), 1, 2);
		switchSubtree(quad, get<0>(c), 2, 1);
		res += "," + printOptimalSubtreeWithSupport(quad, get<1>(c), get<0>(c), support, lambda, qInfo, weight);
		if (support == 3) {
			ofstream fcsv("freqQuad.csv");
			printFreqQuadCSV(get<0>(c), get<1>(c), names[0], fcsv, qInfo);
			printFreqQuadCSV(get<1>(c), get<0>(c), names[0], fcsv, qInfo);
		}
		return res + ")," + names[0] + ");";
	}
#endif
};

#ifdef SUPPORT
const string HELP_TEXT_1 = "binary_path [-c constraintSubtreeFilePath -g guideTreeFilePath -o oFilePath -r nRound -s nSample -p probability -t nThread -u supportLevel";
const string HELP_TEXT_2 = R"V0G0N(] inputList
-c  path to constraint subtree file
-g  path to guide tree file
-o  path to output file (default: stdout)
-r  number of total rounds of placements (default: 5)
-s  number of total rounds of subsampling (default: 0)
-p  subsampling probability of keeping each taxon (default: 0.5)
-t  number of threads (default: 1)
-l  rate lambda of Yule process under which the species tree is modeled (default: 0.5)
-u  output support level (0: no output support value, 1<default>: branch local posterior probability, 2: detailed, 3: freqQuad.csv)
)V0G0N";
#else
const string HELP_TEXT_1 = "binary_path [-c constraintSubtreeFilePath -g guideTreeFilePath -o oFilePath -r nRound -s nSample -p probability -t nThread";
const string HELP_TEXT_2 = R"V0G0N(] inputList
-c  path to constraint subtree file
-g  path to guide tree file
-o  path to output file (default: stdout)
-r  number of total rounds of placements (default: 5)
-s  number of total rounds of subsampling (default: 0)
-p  subsampling probability of keeping each taxon (default: 0.5)
-t  number of threads (default: 1)
)V0G0N";
#endif

struct MetaAlgorithm{
	struct Option{
		unordered_set<string> options;
		bool check(const char* c1, const char* c2){
			options.insert(c2);
			return strcmp(c1, c2) == 0;
		}
		bool isValid(char* c2){
			return options.count(c2);
		}
	};

	vector<string> files, names;
	int nThreads = 1, nRounds = 4, nSample = 4, nBatch = 8, fold = 0, nThread1, nThread2 = 1, support = 1;
	double p = 0.5, lambda = 0.5;
	string outputFile, guideFile, constraintFile, constraintTree;
	ofstream fileOut;
	unordered_map<string, int> name2id;
	Option opt;
	
	TripartitionInitializer tripInit;
	vector<TripartitionInitializer> batchInit;
	
	MetaAlgorithm(){}
	
	MetaAlgorithm(int argc, char** argv, string helpTextS1, string helpTextS2){
		initialize(argc, argv, helpTextS1, helpTextS2);
	}
	
	void initialize(int argc, char** argv, string helpTextS1, string helpTextS2){
		string version = ALG_VERSION;
		#ifdef OBJECTIVE_VERSION
			version = version + "." + OBJECTIVE_VERSION;
		#endif
		#ifdef DRIVER_VERSION
			version = version + "." + DRIVER_VERSION;
		#endif
		#ifdef TUTORIAL
			MDGenerator::version = version;
		#endif
		
	#ifdef ARG_PARSER
		cerr << ARG.getFullName() << endl;
		cerr << "Version: " << version << endl;
		ARG.addStringArg('c', "constraint", "", "Newick file containing a binary species tree to place missing species on");
		ARG.addStringArg('g', "guide", "", "Newick file containing binary trees as guide trees");
		ARG.addStringArg('o', "output", "<standard output>", "File name for the output species tree", true);
		ARG.addIntArg('r', "round", 4, "Number of initial rounds of placements");
		ARG.addIntArg('s', "subsample", 4, "Number of rounds of subsampling per exploration step");
		ARG.addDoubleArg('p', "proportion", 0.5, "Proportion of taxa in the subsample in naive algorithm");
		ARG.addIntArg('t', "thread", 1, "Number of threads", true);
		ARG.addIntArg(0, "seed", 233, "Seed for pseudorandomness");
		ARG.addFlag('C', "scoring", "Scoring the full species tree file after `-c` without exploring other topologies (`-r 1 -s 0`)", [&](){
			ARG.getIntArg("round") = 1; ARG.getIntArg("subsample") = 0;
		}, true);
		ARG.addFlag('R', "moreround", "More rounds of placements and subsampling (`-r 16 -s 16`)", [&](){
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
		
		ARG.parse(argc, argv);
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
	#else
		cerr << "Version: " << version << endl;
		if (argc == 1) {cerr << HELP_TEXT_1 << helpTextS1 << HELP_TEXT_2 << helpTextS2; exit(0);}
		for (int i = 1; i < argc; i += 2){
			if (opt.check(argv[i], "-y")) {i--; continue;}
			
			if (opt.check(argv[i], "-c")) constraintFile = argv[i + 1];
			if (opt.check(argv[i], "-g")) guideFile = argv[i + 1];
			if (opt.check(argv[i], "-o")) outputFile = argv[i + 1];
			if (opt.check(argv[i], "-r")) sscanf(argv[i + 1], "%d", &nRounds);
			if (opt.check(argv[i], "-s")) sscanf(argv[i + 1], "%d", &nSample);
			if (opt.check(argv[i], "-p")) sscanf(argv[i + 1], "%lf", &p);
			if (opt.check(argv[i], "-t")) sscanf(argv[i + 1], "%d", &nThreads);
			if (opt.check(argv[i], "-l")) sscanf(argv[i + 1], "%lf", &lambda);
			if (opt.check(argv[i], "-u")) sscanf(argv[i + 1], "%d", &support);
			if (opt.check(argv[i], "-h")) {cerr << HELP_TEXT_1 << helpTextS1 << HELP_TEXT_2 << helpTextS2; exit(0);}
		}
	#endif
		#ifdef USE_CUDA
		nThread2 = 1;
		nThread1 = 1;
		#else
		nThread2 = nThreads;
		nThread1 = 1;
		#endif
		TP.initialize(nThreads);
		batchInit.resize(nBatch);
	}
	
	void searchSpace(ConstrainedOptimizationAlgorithm &alg, const vector<int> &batchId, const vector<int> &seeds){
		for (int i = 0; i < batchId.size(); i++){
			ConstrainedOptimizationAlgorithm batchAlg(names.size(), batchInit[batchId[i]], names, seeds[i]);
			string tree = batchAlg.run(1, 1).second;
			alg.addGuideTree(tree, name2id);
		}
	}
	
	pair<score_t, string> generateSearchSpace(ConstrainedOptimizationAlgorithm &alg){
		if (nRounds > 0){
			vector<vector<int> > batchId(nThreads), seeds(nThreads);
			vector<thread> thrds;
			
			for (int i = 0; i < nRounds * fold; i++){
				batchId[i % nThreads].push_back(i % nBatch);
				seeds[i % nThreads].push_back(rand());
			}
			
			for (int i = 1; i < nThreads; i++){
				thrds.emplace_back(&MetaAlgorithm::searchSpace, this, ref(alg), ref(batchId[i]), ref(seeds[i]));
			}
			
			searchSpace(alg, batchId[0], seeds[0]);
			
			for (int i = 0; i < nThreads - 1; i++){
				thrds[i].join();
			}
			cerr << "*** End of Batching ***" << endl;
		}
		return alg.run(nRounds, nThread1);
	}
	
	pair<score_t, string> run(){
		ostream &fout = (outputFile == "" || outputFile == "<standard output>") ? cout : fileOut;
		if (outputFile != "" && outputFile != "<standard output>") fileOut.open(outputFile);
		
		cerr << "#Species: " << names.size() << endl;
		cerr << "#Rounds: " << nRounds << endl;
		cerr << "#Samples: " << nSample << endl;
		cerr << "#Threads: " << nThread1 << "x" << nThread2 << endl;
		cerr << "p = " << p << endl;
		
		ConstrainedOptimizationAlgorithm alg(names.size(), tripInit, names);
		
		if (guideFile != ""){
			ifstream fin(guideFile);
			string tree;
			while (getline(fin, tree)){
				if (tree.size()) alg.addGuideTree(tree, name2id);
			}
		}
		
		if (constraintFile != ""){
			ifstream fin(constraintFile);
			getline(fin, constraintTree);
		}
		else{
			ROUND_NN = 20 + 2 * sqrt(names.size()) * log2(names.size());
		}
		
		auto res = (constraintTree == "") ? generateSearchSpace(alg) : alg.constrainedRun(nRounds, nThreads, constraintTree, name2id);
		cerr << "Initial score: " << (double) res.first << endl;
		cerr << "Initial tree: " << res.second << endl;
		if (constraintTree == "") {
			cerr << "*** Subsample Process ***" << endl;
			score_t prevS;
			int roundNum = 0;
			do {
				prevS = res.first;
				res = alg.run(nSample, nThread1, p);
				cerr << "Current score: " << (double) res.first << endl;
				cerr << "Current tree: " << res.second << endl;
				roundNum++;
			}
			while (prevS + ERROR_TOLERANCE < res.first && roundNum < 20);
			if (prevS + ERROR_TOLERANCE < res.first) cerr << "Search stopped due to excessive rounds. Tree may not be optimal!" << endl;
		}
		
		string output = res.second;
		cerr << "Final Tree: " << output << endl;
		#if defined(SUPPORT) || defined(G_SUPPORT)
		double w = 1;
		#ifdef SUPPORT
		w = ARG.getDoubleArg("downweightrepeat");
		#endif
		if (support) output = alg.printOptimalTreeWithSupport(support, lambda, w);
		#endif
		fout << output << endl;
		
		return res;
	}
};

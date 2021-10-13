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

#ifdef SUPPORT
#include "incbeta.c"
#include <cmath>
#endif

using namespace std;

struct Hasher{
	size_t operator()(const hash_t h) const{
		return (size_t) h;
	}
};

const int ROUND_NN = 50;

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
	
	void switchSubtree(int v, int src, int tgt){
		if (leafId(v) != -1){
			//trip.rmv(src, leafId(v));
			//trip.add(tgt, leafId(v));
			trip.update(tgt, leafId(v));
		}
		else{
			switchSubtree(light(v), src, tgt);
			switchSubtree(heavy(v), src, tgt);
		}
	}
	
	void scoringPlacementDP(int v, int i){
		if (leafId(v) != -1){
			//Wv|-|i
			sNo(v) = 0;
			//W|v|i
			trip.update(1, leafId(v));
			sYes(v) = trip.score();
			placeHeavy(v) = false;
			placeLight(v) = false;
		}
		else{
			//AC|B|i
			scoringPlacementDP(light(v), i);
			//ABC|-|i
			switchSubtree(light(v), 1, 0);
			//BC|A|i
			scoringPlacementDP(heavy(v), i);
			//C|A|Bi
			switchSubtree(light(v), 0, 2);
			score_t up, left, right, sibling;
			right = trip.score();
			//C|Ai|B
			trip.update(1, i);
			left = trip.score();
			//Ci|A|B
			trip.update(0, i);
			up = trip.score();
			//C|A|Bi
			trip.update(2, i);
			//C|AB|i
			switchSubtree(light(v), 2, 1);
			sibling = trip.score();
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
	
	int locateBranchAdding(int i){
		//trip.addTotal(i);
		trip.update(2, i);
		trip.update(0, rootLeafId);
		switchSubtree(rootNodeId, -1, 0);
		scoringPlacementDP(rootNodeId, i);
		int v = rootNodeId;
		while (placeHeavy(v) || placeLight(v)){
			if (placeHeavy(v)) v = heavy(v);
			else v = light(v);
		}
		return v;
	}
	
	void place(int i){
		addSibling(locateBranchAdding(i), i);
	}
	
	int locateBranch(int i){
		int v = locateBranchAdding(i);
		trip.update(-1, i); //trip.rmvTotal(i);
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
			switchSubtree(heavy(v), 1, 0);
			switchSubtree(c, 2, 1);
			abe_c_d = trip.score();
			//ABE|CD|-
			switchSubtree(d, 2, 1);
			if (recurse) abe_c_d = nnLocal(light(v), false, abe_c_d);
			//E|AB|CD
			switchSubtree(heavy(v), 0, 1);
			switchSubtree(light(v), 1, 2);
		}
		if (leafId(heavy(v)) == -1){
			int b = light(heavy(v));
			//CDE|A|B
			switchSubtree(light(v), 2, 0);
			switchSubtree(b, 1, 2);
			cde_a_b = trip.score();
			//CDE|AB|-
			switchSubtree(b, 2, 1);
			if (recurse) cde_a_b = nnLocal(heavy(v), false, cde_a_b);
			//E|AB|CD
			switchSubtree(light(v), 0, 2);
		}
		int best_case = 0;
		score_t best_score = abe_c_d + cde_a_b + e_ab_cd, top_score = e_ab_cd;
		if (leafId(light(v)) == -1){
			int c = heavy(light(v)), d = light(light(v));
			//DE|AB|C
			switchSubtree(d, 2, 0);
			score_t de_ab_c = trip.score();
			//E|ABD|C
			switchSubtree(d, 0, 1);
			score_t e_abd_c = trip.score();
			//E|ABC|D
			switchSubtree(d, 1, 2);
			switchSubtree(c, 2, 1);
			score_t e_abc_d = trip.score();
			//CE|AB|D
			switchSubtree(c, 1, 0);
			score_t ce_ab_d = trip.score();
			
			if (best_score < cde_a_b + de_ab_c + e_abc_d) {
				best_score = cde_a_b + de_ab_c + e_abc_d;
				best_case = 1;
				top_score = e_abc_d;
			}
			if (best_score < cde_a_b + ce_ab_d + e_abd_c) {
				best_score = cde_a_b + ce_ab_d + e_abd_c;
				best_case = 2;
				top_score = e_abd_c;
			}
			//E|AB|CD
			switchSubtree(c, 0, 2);
		}
		if (leafId(heavy(v)) == -1){
			int a = heavy(heavy(v)), b = light(heavy(v));
			//BE|A|CD
			switchSubtree(b, 1, 0);
			score_t be_a_cd = trip.score();
			//E|A|BCD
			switchSubtree(b, 0, 2);
			score_t e_a_bcd = trip.score();
			//E|ACD|B
			switchSubtree(light(v), 2, 1);
			score_t e_acd_b = trip.score();
			//AE|B|CD
			switchSubtree(light(v), 1, 2);
			switchSubtree(a, 1, 0);
			switchSubtree(b, 2, 1);
			score_t ae_b_cd = trip.score();
			
			if (best_score < abe_c_d + ae_b_cd + e_a_bcd) {
				best_score = abe_c_d + ae_b_cd + e_a_bcd;
				best_case = 3;
				top_score = e_a_bcd;
			}
			if (best_score < abe_c_d + be_a_cd + e_acd_b) {
				best_score = abe_c_d + be_a_cd + e_acd_b;
				best_case = 4;
				top_score = e_acd_b;
			}
			//E|AB|CD
			switchSubtree(a, 0, 1);
		}
		//E|ABCD|-
		switchSubtree(light(v), 2, 1);
		return nnPerform(v, best_case, top_score);
	}
	
	score_t nnMove(int v){
		if (leafId(v) != -1){
			trip.update(1, leafId(v));
			return 0;
		}
		else{
			//ABE|CD|-
			score_t abe_c_d = nnMove(light(v));
			score_t ae_b_cd = 0;
			if (leafId(heavy(v)) == -1){
				int b = light(heavy(v));
				//AE|B|CD
				switchSubtree(b, 0, 1);
				switchSubtree(light(v), 1, 2);
				ae_b_cd = trip.score();
				//ABCDE|-|-
				switchSubtree(b, 1, 0);
				switchSubtree(light(v), 2, 0);
			}
			//ABCDE|-|-
			else switchSubtree(light(v), 1, 0);
			//CDE|AB|-
			score_t cde_a_b = nnMove(heavy(v));
			//E|AB|CD
			switchSubtree(light(v), 0, 2);
			score_t e_ab_cd = trip.score();
			int best_case = 0;
			score_t best_score = abe_c_d + cde_a_b + e_ab_cd, top_score = e_ab_cd;
			
			if (leafId(light(v)) == -1){
				int c = heavy(light(v)), d = light(light(v));
				//DE|AB|C
				switchSubtree(d, 2, 0);
				score_t de_ab_c = trip.score();
				//E|ABD|C
				switchSubtree(d, 0, 1);
				score_t e_abd_c = trip.score();
				//E|ABC|D
				switchSubtree(d, 1, 2);
				switchSubtree(c, 2, 1);
				score_t e_abc_d = trip.score();
				//CE|AB|D
				switchSubtree(c, 1, 0);
				score_t ce_ab_d = trip.score();
				
				if (best_score < cde_a_b + de_ab_c + e_abc_d) {
					best_score = cde_a_b + de_ab_c + e_abc_d;
					best_case = 1;
					top_score = e_abc_d;
				}
				if (best_score < cde_a_b + ce_ab_d + e_abd_c) {
					best_score = cde_a_b + ce_ab_d + e_abd_c;
					best_case = 2;
					top_score = e_abd_c;
				}
				//E|AB|CD
				switchSubtree(c, 0, 2);
			}
			if (leafId(heavy(v)) == -1){
				int b = light(heavy(v));
				//E|A|BCD
				switchSubtree(b, 1, 2);
				score_t e_a_bcd = trip.score();
				//BE|A|CD
				switchSubtree(b, 2, 0);
				score_t be_a_cd = trip.score();
				//E|ACD|B
				switchSubtree(b, 0, 2);
				switchSubtree(light(v), 2, 1);
				score_t e_acd_b = trip.score();
				
				if (best_score < abe_c_d + ae_b_cd + e_a_bcd) {
					best_score = abe_c_d + ae_b_cd + e_a_bcd;
					best_case = 3;
					top_score = e_a_bcd;
				}
				if (best_score < abe_c_d + be_a_cd + e_acd_b) {
					best_score = abe_c_d + be_a_cd + e_acd_b;
					best_case = 4;
					top_score = e_acd_b;
				}
				//E|AB|CD
				switchSubtree(b, 2, 1);
				switchSubtree(light(v), 1, 2);
			}
			//E|ABCD|-
			switchSubtree(light(v), 2, 1);
			return nnPerform(v, best_case, top_score);
		}
	}
	
	pair<bool, hash_t> tripHashGenerator(int v){
		if (leafId(v) != -1){
			//trip.reset();
			//trip.rmv(0, leafId(v));
			//trip.add(1, leafId(v));
			trip.update(1, leafId(v));
			return {leafId(v) == 0, taxonHash[leafId(v)]};
		}
		else{
			//AC|B|-
			pair<bool, hash_t> hLight = tripHashGenerator(light(v));
			//ABC|-|-
			switchSubtree(light(v), 1, 0);
			//BC|A|-
			pair<bool, hash_t> hHeavy = tripHashGenerator(heavy(v));
			//C|A|B
			switchSubtree(light(v), 0, 2);
			score_t s = trip.score();
			//C|AB|-
			switchSubtree(light(v), 2, 1);
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
		if (rNN >= 0) {
			trip.update(0, rootLeafId);
			switchSubtree(rootNodeId, -1, 0);
			nnMove(rootNodeId);
			cerr << "#localmove:" << ROUND_NN - rNN << "/" << ROUND_NN << endl;
		}
		trip.update(0, rootLeafId);
		switchSubtree(rootNodeId, -1, 0);
		tripHashGenerator(rootNodeId);
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
	
	ConstrainedOptimizationAlgorithm(const int ntaxa, const TripartitionInitializer &tripInit, const vector<string> &names, const int seed = 2333): ntaxa(ntaxa), tripInit(tripInit), names(names), generator(seed){
		taxonHash.push_back(0);
		for (int i = 1; i < ntaxa; i++){
			hash_t r = randomHash(generator);
			taxonHash.push_back(r);
			hash[r] = nodes.size();
			nodes.emplace_back(i, r);
			taxonHash[0] -= r;
		}
	}
	
	ConstrainedOptimizationAlgorithm(const ConstrainedOptimizationAlgorithm &alg): ntaxa(alg.ntaxa), tripInit(alg.tripInit), names(alg.names), nodes(alg.nodes), hash(alg.hash), taxonHash(alg.taxonHash){}
	
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
	
	int guildSubtree(PlacementAlgorithm &pAlg, const string &tree, const unordered_map<string, int> &name2id, int &i, unordered_set<int> &added){
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
				int nodeL = guildSubtree(pAlg, tree, name2id, i, added);
				int nodeR = guildSubtree(pAlg, tree, name2id, i, added);
				if (nodeL != -1){
					if (nodeR != -1) ret = pAlg.makeNode(nodeL, nodeR, false);
					else ret = nodeL;
				}
				else ret = nodeR;
			}
			else {
				int nodeL = guildSubtree(pAlg, tree, name2id, i, added);
				if (pAlg.rootLeafId == -1) ret = guildSubtree(pAlg, tree, name2id, i, added);
				else if (nodeL == -1) ret = pAlg.rootNodeId = guildSubtree(pAlg, tree, name2id, i, added);
				else {
					int nodeR = guildSubtree(pAlg, tree, name2id, i, added);
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
		guildSubtree(pAlg, tree, name2id, i, added);
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
			jobs.push_back(createPlacementAlgorithm(subsampleRate));
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
			//pAlg.trip.rmvTotal(i);
			pAlg.trip.update(-1, i);
			return false;
		}
	}
	
	void twoStepWorkflow(PlacementAlgorithm &pAlg){
		int N, rId;
		{ const lock_guard<mutex> lock(mtx); N = ntaxa; rId = roundId; }
		int n = sqrt(N);
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
				for (int i = 0; i < n; i++) pAlg.trip.update(-1, order[i]); //pAlg.trip.rmvTotal(order[i]);
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
			cerr << "Subsampled nodes: " << pAlg.nodes.size() << endl;
		}
		{
			vector<int> abnormalOrder;
			vector<PlacementAlgorithm::Node> pNodes(pAlg.nodes);
			int rootNodeId = pAlg.rootNodeId;//, rootLeafId = pAlg.rootLeafId;
			{
				vector<PlacementAlgorithm::Node> originalNodes(pAlg.nodes);
				int originalRootNodeId = pAlg.rootNodeId;//, originalRootLeafId = pAlg.rootLeafId;
				
				vector<vector<int> > nodeBranch(pNodes.size());
				cerr << "pNodes: " << pNodes.size() << endl;
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
					for (int i: added) pAlg.trip.update(-1, i); //pAlg.trip.rmvTotal(i);
				}
			}
			//for (int i = n; i < N; i++) pAlg.trip.addTotal(order[i]);
			//for (int i: abnormalOrder) pAlg.trip.rmvTotal(i);
			pAlg.order = abnormalOrder;
			pAlg.nodes = pNodes;
			pAlg.rootNodeId = rootNodeId;
			//pAlg.rootLeafId = rootLeafId;
			pAlg.orderId = 0;
		}
		cerr << "Remaining: " << pAlg.order.size() << endl;
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

#ifdef SUPPORT
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
			unordered_map<int, tuple<array<double, 3>, array<double, 3>, string> > &qInfo){
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
		res += printOptimalSubtreeWithSupport(quad, get<0>(c), get<1>(c), support, lambda, qInfo) + ",";
		//ru|c0|c1|-
		switchSubtree(quad, get<0>(c), 2, 1);
		switchSubtree(quad, get<1>(c), 1, 2);
		res += printOptimalSubtreeWithSupport(quad, get<1>(c), get<0>(c), support, lambda, qInfo) + ")";
		//r|u|c1|c0
		switchSubtree(quad, u, 0, 1);
		switchSubtree(quad, get<0>(c), 1, 3);
		array<double, 3> score = quad.score();
		//r|u|c0c1|-
		switchSubtree(quad, get<0>(c), 3, 2);
		if (score[0] + score[1] + score[2] < 1e-8) {
			if (support == 1) return res + to_string(1.0 / 3.0);
			else return res + "'support=(0,0,0);p=(0.333,0.333,0.333)'"; 
		}
		double tscore = score[0] + score[1] + score[2];
		double i0 = 1.0 - incbeta(score[0] + 1.0, tscore + lambda * 2 - score[0], 1.0 / 3.0);
		double i1 = 1.0 - incbeta(score[1] + 1.0, tscore + lambda * 2 - score[1], 1.0 / 3.0);
		double i2 = 1.0 - incbeta(score[2] + 1.0, tscore + lambda * 2 - score[2], 1.0 / 3.0);
		double lb0 = lgamma(score[0] + 1.0) + lgamma(tscore - score[0] + lambda * 2) - lgamma(tscore + 1.0 + lambda * 2);
		double lb1 = lgamma(score[1] + 1.0) + lgamma(tscore - score[1] + lambda * 2) - lgamma(tscore + 1.0 + lambda * 2);
		double lb2 = lgamma(score[2] + 1.0) + lgamma(tscore - score[2] + lambda * 2) - lgamma(tscore + 1.0 + lambda * 2);
		if (support == 1) res += to_string(i0 / (i0 + i1 * exp(log(2.0) * (score[1] - score[0]) + lb1 - lb0) + i2 * exp(log(2.0) * (score[2] - score[0]) + lb2 - lb0)));
		else {
			array<double, 3> p;
			p[0] = i0 / (i0 + i1 * exp(log(2.0) * (score[1] - score[0]) + lb1 - lb0) + i2 * exp(log(2.0) * (score[2] - score[0]) + lb2 - lb0));
			p[1] = i1 / (i1 + i0 * exp(log(2.0) * (score[0] - score[1]) + lb0 - lb1) + i2 * exp(log(2.0) * (score[2] - score[1]) + lb2 - lb1));
			p[2] = i2 / (i2 + i1 * exp(log(2.0) * (score[1] - score[2]) + lb1 - lb2) + i0 * exp(log(2.0) * (score[0] - score[2]) + lb0 - lb2));
			res += "'support=(" + to_string(score[0]) + "," + to_string(score[1]) + "," + to_string(score[2]) + ");p=(";
			res += to_string(p[0]) + "," + to_string(p[1]) + "," + to_string(p[2]) + ")'";
			if (support == 3) qInfo[v] = make_tuple(score, p, get<2>(qInfo[get<0>(c)]) + "," + get<2>(qInfo[get<1>(c)]));
		}
		if (3 * score[0] > tscore) res += ":" + to_string(max(0.0, -log(1.5 - 1.5 * score[0] / (tscore + lambda * 2))));
		else res += ":0";
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
	
	string printOptimalTreeWithSupport(int support, double lambda){
		unordered_map<int, tuple<array<double, 3>, array<double, 3>, string> > qInfo;
		string res;
		Quadrupartition quad(tripInit);
		quad.update(0, 0);
		for (int i = 1; i < names.size(); i++) quad.update(2, i);
		int v = hash[-taxonHash[0]];
		tuple<int, int, score_t> c = nodes[v].children[nodes[v].bestChild];
		//0|c1|c0|-
		switchSubtree(quad, get<1>(c), 2, 1);
		res = "((" + printOptimalSubtreeWithSupport(quad, get<0>(c), get<1>(c), support, lambda, qInfo);
		//0|c0|c1|-
		switchSubtree(quad, get<1>(c), 1, 2);
		switchSubtree(quad, get<0>(c), 2, 1);
		res += "," + printOptimalSubtreeWithSupport(quad, get<1>(c), get<0>(c), support, lambda, qInfo);
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
const string HELP_TEXT_1 = "feast [-c constraintSubtreeFilePath -g guideTreeFilePath -o oFilePath -r nRound -s nSample -p probability -t nThread -u supportLevel";
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
	vector<string> files, names;
	int nThreads = 1, nRounds = 4, nSample = 4, nBatch = 8, fold = 0, nThread1, nThread2 = 1, support = 1;
	double p = 0.5, lambda = 0.5;
	string outputFile, guideFile, constraintFile, constraintTree;
	ofstream fileOut;
	unordered_map<string, int> name2id;
	
	TripartitionInitializer tripInit;
	vector<TripartitionInitializer> batchInit;
	
	MetaAlgorithm(){}
	
	MetaAlgorithm(int argc, char** argv, string helpTextS1, string helpTextS2){
		initialize(argc, argv, helpTextS1, helpTextS2);
	}
	
	void initialize(int argc, char** argv, string helpTextS1, string helpTextS2){
		if (argc == 1) {cerr << HELP_TEXT_1 << helpTextS1 << HELP_TEXT_2 << helpTextS2; exit(0);}
		for (int i = 1; i < argc; i += 2){
			if (strcmp(argv[i], "-y") == 0) {i--; continue;}
			
			if (strcmp(argv[i], "-c") == 0) constraintFile = argv[i + 1];
			if (strcmp(argv[i], "-g") == 0) guideFile = argv[i + 1];
			if (strcmp(argv[i], "-o") == 0) outputFile = argv[i + 1];
			if (strcmp(argv[i], "-r") == 0) sscanf(argv[i + 1], "%d", &nRounds);
			if (strcmp(argv[i], "-s") == 0) sscanf(argv[i + 1], "%d", &nSample);
			if (strcmp(argv[i], "-p") == 0) sscanf(argv[i + 1], "%lf", &p);
			if (strcmp(argv[i], "-t") == 0) sscanf(argv[i + 1], "%d", &nThreads);
			if (strcmp(argv[i], "-l") == 0) sscanf(argv[i + 1], "%lf", &lambda);
			if (strcmp(argv[i], "-u") == 0) sscanf(argv[i + 1], "%d", &support);
			if (strcmp(argv[i], "-h") == 0) {cerr << HELP_TEXT_1 << helpTextS1 << HELP_TEXT_2 << helpTextS2; exit(0);}
		}
		
		if (nRounds < nThreads && nRounds > 0){
			nThread2 = nThreads / nRounds;
			nThread1 = nRounds;
		}
		else {
			nThread2 = 1;
			nThread1 = nThreads;
		}
		
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
		ostream &fout = (outputFile == "") ? cout : fileOut;
		if (outputFile != "") fileOut.open(outputFile);
		
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
				alg.addGuideTree(tree, name2id);
			}
		}
		
		if (constraintFile != ""){
			ifstream fin(constraintFile);
			getline(fin, constraintTree);
		}
		
		auto res = (constraintTree == "") ? generateSearchSpace(alg) : alg.constrainedRun(nRounds, nThreads, constraintTree, name2id);
		if (constraintTree == "") {
			cerr << "*** Subsample Process ***" << endl;
			score_t prevS;
			int roundNum = 0;
			do {
				prevS = res.first;
				res = alg.run(nSample, nThread1, p);
				cerr << "Current score: " << (double) res.first << endl;
				roundNum++;
			}
			while (prevS < res.first && roundNum < 5);
		}
		
		string output = res.second;
		cerr << "Final Tree: " << output << endl;
		#ifdef SUPPORT
		if (support) output = alg.printOptimalTreeWithSupport(support, lambda);
		#endif
		fout << output << endl;
		
		return res;
	}
};

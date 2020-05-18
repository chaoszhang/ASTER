#include "sequence.hpp"

typedef __int64 hash_t;

using namespace std;

struct PlacementAlgorithm{
	struct Node{
		Node(int leafId): leafId(leafId){}
		
		int parent = -1, heavy = -1, light = -1; //parent 0, heavy 1, light 2
		int leafId = -1, leafCnt = 1;
		score_t sYes = 0, sNo = 0;
		bool placeHeavy = false, placeLight = false; 
	};
	
	vector<hash_t> taxonHash; // sum of hash = 0
	vector<int> order;
	vector<Node> nodes;
	int rootNodeId = 0, rootLeafId = 0, orderId = 0;
	
	Tripartition trip;
	
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
	
	void addSibling(int v, int i){
		int u = nodes.size();
		nodes.emplace_back(i);
		int w = nodes.size();
		nodes.emplace_back(-1);
		parent(w) = parent(v);
		if (parent(v) == -1) rootNodeId = w;
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
	}
	
	void defaultInitializer(){
		orderId = 3;
		nodes.emplace_back(order[1]);
		rootNodeId = 0;
		rootLeafId = order[0];
		trip.addTotal(order[0]);
		trip.addTotal(order[1]);
		addSibling(0, order[2]);
	}
	
	void switchSubtree(int v, int src, int tgt){
		if (leafId(v) != -1){
			trip.rmv(src, leafId(v));
			trip.add(tgt, leafId(v));
		}
		else{
			switchSubtree(light(v), src, tgt);
			switchSubtree(heavy(v), src, tgt);
		}
	}
	
	void scoringPlacementDP(int v, int i){
		if (leafId(v) != -1){
			trip.reset();
			trip.rmv(0, leafId(v));
			trip.add(1, leafId(v));
			sNo(v) = 0;
			trip.rmv(0, i);
			trip.add(2, i);
			sYes(v) = trip.score();
			placeHeavy(v) = false;
			placeLight(v) = false;
			trip.rmv(2, i);
			trip.add(0, i);
		}
		else{
			//BCi|A|-
			scoringPlacementDP(light(v), i);
			scoringPlacementDP(heavy(v), i);
			//Ci|A|B
			switchSubtree(light(v), 0, 2);
			score_t up, left, right, sibling;
			up = trip.score();
			//C|Ai|B
			trip.rmv(0, i);
			trip.add(1, i);
			left = trip.score();
			//C|A|Bi
			trip.rmv(1, i);
			trip.add(2, i);
			right = trip.score();
			//C|AB|i
			switchSubtree(light(v), 2, 1);
			sibling = trip.score();
			//Ci|AB|-
			trip.rmv(2, i);
			trip.add(0, i);
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
				placeHeavy(v) = false;
				placeLight(v) = true;
			}
		}
	}
	
	void place(int i){
		trip.addTotal(i);
		scoringPlacementDP(rootNodeId, i);
		int v = rootNodeId;
		while (placeHeavy(v) || placeLight(v)){
			if (placeHeavy(v)) v = heavy(v);
			else v = light(v);
		}
		addSibling(v, i);
	}
	
	void run(){
		defaultInitializer();
		for (; orderId < taxonHash.size(); orderId++){
			place(order[orderId]);
		}
	}
};

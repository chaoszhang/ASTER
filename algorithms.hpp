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

using namespace std;

struct Hasher{
	size_t operator()(const hash_t h) const{
		return (size_t) h;
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
	int rootNodeId = -1, rootLeafId = -1, orderId = 0;
	
	Tripartition trip;
	
	PlacementAlgorithm(const vector<hash_t> &taxonHash, const TripartitionInitializer &tripInit): taxonHash(taxonHash), trip(tripInit){}
	
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
	
	void defaultInitializer(){
		if (rootLeafId == -1){
			rootLeafId = order[orderId];
			trip.addTotal(order[orderId]);
			orderId++;
		}
		if (rootNodeId == -1){
			rootNodeId = nodes.size();
			nodes.emplace_back(order[orderId]);
			trip.addTotal(order[orderId]);
			orderId++;
		}
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
			//Wvi|-|-
			trip.reset();
			//Wi|v|-
			trip.rmv(0, leafId(v));
			trip.add(1, leafId(v));
			sNo(v) = 0;
			//W|v|i
			trip.rmv(0, i);
			trip.add(2, i);
			sYes(v) = trip.score();
			placeHeavy(v) = false;
			placeLight(v) = false;
			//Wi|v|-
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
				placeHeavy(v) = true;
				placeLight(v) = false;
			}
		}
	}
	
	int locateBranchAdding(int i){
		trip.addTotal(i);
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
		trip.rmvTotal(i);
		return v;
	}
	
	pair<bool, hash_t> tripHashGenerator(int v){
		if (leafId(v) != -1){
			trip.reset();
			trip.rmv(0, leafId(v));
			trip.add(1, leafId(v));
			return {leafId(v) == 0, taxonHash[leafId(v)]};
		}
		else{
			pair<bool, hash_t> hLight = tripHashGenerator(light(v)), hHeavy = tripHashGenerator(heavy(v));
			//BC|A|- -> C|A|B
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
			cerr << "Placing " << orderId << "/" << order.size() << endl;
			place(order[orderId]);
		}
		tripHashGenerator(rootNodeId);
	}
	
	
	string printSubtree(int v){
		if (leafId(v) != -1) return to_string(leafId(v));
		return string("(") + printSubtree(heavy(v)) + "," + printSubtree(light(v)) + ")";
	}
	
	string printTree(){
		return string("(") + printSubtree(rootNodeId) + "," + to_string(rootLeafId) + ");";
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
				pAlg.trip.addTotal(nodes[v].leafId);
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
				pAlg.trip.addTotal(0);
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
				pAlg.trip.addTotal(nodes[v].leafId);
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
				pAlg.trip.addTotal(0);
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
	
	int guildSubtree(PlacementAlgorithm &pAlg, const string &tree, const unordered_map<string, int> &name2id, int &i){
		int ret;
		if (tree[i] != '('){
			string s;
			for (; tree[i] != ':' && tree[i] != ')' && tree[i] != ',' && tree[i] != ';'; i++){
				if (tree[i] != '\"' && tree[i] != '\'') s += tree[i];
			}
			int id = name2id.at(s);
			pAlg.trip.addTotal(id);
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
				int nodeL = guildSubtree(pAlg, tree, name2id, i);
				int nodeR = guildSubtree(pAlg, tree, name2id, i);
				if (nodeL != -1){
					if (nodeR != -1) ret = pAlg.makeNode(nodeL, nodeR, false);
					else ret = nodeL;
				}
				else ret = nodeR;
			}
			else {
				int nodeL = guildSubtree(pAlg, tree, name2id, i);
				if (pAlg.rootLeafId == -1) ret = guildSubtree(pAlg, tree, name2id, i);
				else if (nodeL == -1) ret = pAlg.rootNodeId = guildSubtree(pAlg, tree, name2id, i);
				else {
					int nodeR = guildSubtree(pAlg, tree, name2id, i);
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
	
	void addGuideTree(const string tree, const unordered_map<string, int> &name2id){
		PlacementAlgorithm pAlg(taxonHash, tripInit);
		int i = 0;
		guildSubtree(pAlg, tree, name2id, i);
		pAlg.run();
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
	
	static void batchWork(vector<PlacementAlgorithm> &jobs, int start, int end){
		for (int i = start; i < end; i++){
			jobs[i].run();
		}
	}
	
	pair<score_t, string> run(int nJobs = 1, int nThrds = 1, double subsampleRate = 0, bool allowTwoStepRun = true){
		if (allowTwoStepRun && ntaxa >= 100) return twoStepRun(nJobs, nThrds);
		vector<PlacementAlgorithm> jobs;
		vector<thread> thrds;
		for (int i = 0; i < nJobs; i++){
			jobs.push_back(createPlacementAlgorithm(subsampleRate));
		}
		for (int i = 1; i < nThrds; i++){
			thrds.emplace_back(batchWork, ref(jobs), i * nJobs / nThrds, (i + 1) * nJobs / nThrds);
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
			pAlg.trip.rmvTotal(i);
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
				for (int i = 0; i < n; i++) pAlg.trip.rmvTotal(order[i]);
				pAlg.tripHash.clear();
				pAlg.order.clear();
				pAlg.nodes.clear();
				pAlg.rootNodeId = -1;
				pAlg.rootLeafId = -1;
				pAlg.orderId = 0;
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
			cerr << "subsampled nodes: " << pAlg.nodes.size() << endl;
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
					cerr << "Binning " << i - n << "/" << N - n << endl; 
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
						cerr << "Placement " << cnt++ << "/" << N - n << endl;
						if (placeNode(pAlg, pNodes, rootNodeId, i, nodeRemap)) added.push_back(i);
						else abnormalOrder.push_back(i);
					}
					for (int i: added) pAlg.trip.rmvTotal(i);
				}
			}
			for (int i = n; i < N; i++) pAlg.trip.addTotal(order[i]);
			for (int i: abnormalOrder) pAlg.trip.rmvTotal(i);
			pAlg.order = abnormalOrder;
			pAlg.nodes = pNodes;
			pAlg.rootNodeId = rootNodeId;
			//pAlg.rootLeafId = rootLeafId;
			pAlg.orderId = 0;
		}
		cerr << "Remaining: " << pAlg.order.size() << endl;
		pAlg.run();
	}
	
	void twoStepBatchWork(vector<PlacementAlgorithm> &jobs, int start, int end){
		for (int i = start; i < end; i++){
			twoStepWorkflow(jobs[i]);
		}
		cerr << "DONE!!!" << endl;
	}
	
	pair<score_t, string> twoStepRun(int nJobs, int nThrds){
		cerr << "Use two-step algorithm!\n";
		//if (roundId != 0) return {computeOptimalTree(), printOptimalTree()};
		vector<PlacementAlgorithm> jobs;
		vector<thread> thrds;
		for (int i = 0; i < nJobs; i++){
			jobs.emplace_back(taxonHash, tripInit);
		}
		for (int i = 1; i < nThrds; i++){
			thrds.emplace_back(&ConstrainedOptimizationAlgorithm::twoStepBatchWork, this, ref(jobs), i * nJobs / nThrds, (i + 1) * nJobs / nThrds);
		}
		twoStepBatchWork(jobs, 0, nJobs / nThrds);
		for (int i = 1; i < nThrds; i++){
			thrds[i - 1].join();
		}
		for (int i = 0; i < nJobs; i++){
			addTripartitions(jobs[i].tripHash);
		}
		
		return {computeOptimalTree(), printOptimalTree()};
	}
};

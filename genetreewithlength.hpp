typedef double score_t;

#include<vector>
#include<array>

using namespace std;

struct TripartitionInitializer{
	struct Node{
		int up = -1, small = -1, large = -1;
		score_t weight = 1;
	};
	
	vector<Node> nodes;
	vector<vector<int> > leafParent;
};

struct Tripartition{
	struct Node{
		score_t x = 0, y = 0, z = 0, tx = 0, ty = 0, tz = 0, q = 0;
		score_t x2a = 0, y2a = 0, z2a = 0, xya = 0, xza = 0, yza = 0;
		score_t x2b = 0, y2b = 0, z2b = 0, xyb = 0, xzb = 0, yzb = 0;
		
		int version = 0, up = -1, small = -1, large = -1; // -1 for dummy!
		score_t weight = 1;
		
		void update(int v, const Node &t){
			if (version != v){
				*this = t;
				version = v;
			}
		}
	};
	
	score_t normal(Node& w){
		Node& u = nodes[w.small];
		Node& v = nodes[w.large];
		u.update(version, totalZ[w.small]);
		v.update(version, totalZ[w.large]);
		
		w.x = (u.x + v.x) * w.weight;
		w.y = (u.y + v.y) * w.weight;
		w.z = (u.z + v.z) * w.weight;
		w.x2a = u.x2a + v.x2a + u.x * v.x;
		w.y2a = u.y2a + v.y2a + u.y * v.y;
		w.z2a = u.z2a + v.z2a + u.z * v.z;
		w.xya = u.xya + v.xya + u.x * v.y + u.y * v.x;
		w.xza = u.xza + v.xza + u.x * v.z + u.z * v.x;
		w.yza = u.yza + v.yza + u.y * v.z + u.z * v.y;
		w.x2b = (u.x2b + v.x2b + u.x * v.x) * w.weight;
		w.y2b = (u.y2b + v.y2b + u.y * v.y) * w.weight;
		w.z2b = (u.z2b + v.z2b + u.z * v.z) * w.weight;
		w.xyb = (u.xyb + v.xyb + u.x * v.y + u.y * v.x) * w.weight;
		w.xzb = (u.xzb + v.xzb + u.x * v.z + u.z * v.x) * w.weight;
		w.yzb = (u.yzb + v.yzb + u.y * v.z + u.z * v.y) * w.weight;
		
		
		w.tx = (u.tx + v.tx + u.y * (v.z2a - v.z2b) + (u.z2a - u.z2b) * v.y + u.z * (v.y2a - v.y2b) + (u.y2a - u.y2b) * v.z) * w.weight;
		w.ty = (u.ty + v.ty + u.x * (v.z2a - v.z2b) + (u.z2a - u.z2b) * v.x + u.z * (v.x2a - v.x2b) + (u.x2a - u.x2b) * v.z) * w.weight;
		w.tz = (u.tz + v.tz + u.x * (v.y2a - v.y2b) + (u.y2a - u.y2b) * v.x + u.y * (v.x2a - v.x2b) + (u.x2a - u.x2b) * v.y) * w.weight;
		score_t oldQ = w.q;
		w.q = u.x * v.tx + u.y * v.ty + u.z * v.tz + v.x * u.tx + v.y * u.ty + v.z * u.tz
			+ u.x2a * v.yza - u.x2b * v.yzb
			+ u.y2a * v.xza - u.y2b * v.xzb
			+ u.z2a * v.xya - u.z2b * v.xyb;
		return w.q - oldQ;
	}
	
	vector<vector<int> > leafParent;
	score_t totalScore = 0;
	vector<Node> nodes, totalZ;
	int version = 0;
	
	Tripartition(TripartitionInitializer init): leafParent(init.leafParent), nodes(init.nodes.size()), totalZ(init.nodes.size()){
		for (int i = 0; i < nodes.size(); i++){
			nodes[i].up = init.nodes[i].up;
			nodes[i].small = init.nodes[i].small;
			nodes[i].large = init.nodes[i].large;
			nodes[i].weight = init.nodes[i].weight;
			totalZ[i] = nodes[i];
		}
	}
	
	void reset(){
		totalScore = 0;
		version++;
	}
	
	void addTotal(int i){
		for (int w: leafParent[i]){
			totalZ[w].z += totalZ[w].weight;
			w = totalZ[w].up;
			while (w != -1){
				Node& u = totalZ[totalZ[w].small];
				Node& v = totalZ[totalZ[w].large];
				totalZ[w].z = (u.z + v.z) * totalZ[w].weight;
				totalZ[w].z2a = u.z2a + v.z2a + u.z * v.z;
				totalZ[w].z2b = (u.z2b + v.z2b + u.z * v.z) * totalZ[w].weight;
				w = totalZ[w].up;
			}
		}
	}
	
	void add(int x, int i){
		for (int u: leafParent[i]){
			nodes[u].update(version, totalZ[u]);
			if (x == 0) nodes[u].z += nodes[u].weight; else if (x == 1) nodes[u].x += nodes[u].weight; else nodes[u].y += nodes[u].weight;
			int w = nodes[u].up;
			while (w != -1){
				nodes[w].update(version, totalZ[w]);
				totalScore += normal(nodes[w]);
				u = w;
				w = nodes[u].up;
			}
		}
	}
	
	void rmv(int x, int i){
		for (int u: leafParent[i]){
			nodes[u].update(version, totalZ[u]);
			if (x == 0) nodes[u].z -= nodes[u].weight; else if (x == 1) nodes[u].x -= nodes[u].weight; else nodes[u].y -= nodes[u].weight;
			int w = nodes[u].up;
			while (w != -1){
				nodes[w].update(version, totalZ[w]);
				totalScore += normal(nodes[w]);
				u = w;
				w = nodes[u].up;
			}
		}
	}
	
	score_t score(){
		return totalScore;
	}
};

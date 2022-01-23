#include<vector>
#include<array>
#include<thread>

using namespace std;

struct TripartitionInitializer{
	struct Node{
		int up = -1, small = -1, large = -1;
		score_t weight = 1;
	};
	
	vector<vector<Node> > nodes;
	vector<vector<vector<int> > > leafParent;
};

struct Tripartition{
	struct Partition{
		struct Node{
			score_t x = 0, y = 0, z = 0, tx = 0, ty = 0, tz = 0, q = 0;
			score_t x2a = 0, y2a = 0, z2a = 0, xya = 0, xza = 0, yza = 0;
			score_t x2b = 0, y2b = 0, z2b = 0, xyb = 0, xzb = 0, yzb = 0;
			
			int up = -1, small = -1, large = -1; // -1 for dummy!
			score_t weight = 1;
			
		};
		
		score_t normal(Node& w){
			Node& u = nodes[w.small];
			Node& v = nodes[w.large];
			//u.update(version, totalZ[w.small]);
			//v.update(version, totalZ[w.large]);
			
			w.x = (u.x + v.x) * 0.5;
			w.y = (u.y + v.y) * 0.5;
			w.z = (u.z + v.z) * 0.5;
			w.x2a = (u.x2a + v.x2a + u.x * v.x) * 0.5;
			w.y2a = (u.y2a + v.y2a + u.y * v.y) * 0.5;
			w.z2a = (u.z2a + v.z2a + u.z * v.z) * 0.5;
			w.xya = (u.xya + v.xya + u.x * v.y + u.y * v.x) * 0.5;
			w.xza = (u.xza + v.xza + u.x * v.z + u.z * v.x) * 0.5;
			w.yza = (u.yza + v.yza + u.y * v.z + u.z * v.y) * 0.5;
			w.x2b = (u.x2b + v.x2b + u.x * v.x) * (1 - w.weight) * 0.5;
			w.y2b = (u.y2b + v.y2b + u.y * v.y) * (1 - w.weight) * 0.5;
			w.z2b = (u.z2b + v.z2b + u.z * v.z) * (1 - w.weight) * 0.5;
			w.xyb = (u.xyb + v.xyb + u.x * v.y + u.y * v.x) * (1 - w.weight) * 0.5;
			w.xzb = (u.xzb + v.xzb + u.x * v.z + u.z * v.x) * (1 - w.weight) * 0.5;
			w.yzb = (u.yzb + v.yzb + u.y * v.z + u.z * v.y) * (1 - w.weight) * 0.5;
			
			
			w.tx = (u.tx + v.tx + u.y * (v.z2a - v.z2b) + (u.z2a - u.z2b) * v.y + u.z * (v.y2a - v.y2b) + (u.y2a - u.y2b) * v.z) * 0.5;
			w.ty = (u.ty + v.ty + u.x * (v.z2a - v.z2b) + (u.z2a - u.z2b) * v.x + u.z * (v.x2a - v.x2b) + (u.x2a - u.x2b) * v.z) * 0.5;
			w.tz = (u.tz + v.tz + u.x * (v.y2a - v.y2b) + (u.y2a - u.y2b) * v.x + u.y * (v.x2a - v.x2b) + (u.x2a - u.x2b) * v.y) * 0.5;
			score_t oldQ = w.q;
			w.q = u.x * v.tx + u.y * v.ty + u.z * v.tz + v.x * u.tx + v.y * u.ty + v.z * u.tz
				+ u.x2a * v.yza - u.x2b * v.yzb
				+ u.y2a * v.xza - u.y2b * v.xzb
				+ u.z2a * v.xya - u.z2b * v.xyb;
			return w.q - oldQ;
		}
		
		vector<vector<int> > leafParent;
		score_t totalScore = 0;
		vector<Node> nodes;
		vector<int> color;
		
		Partition(TripartitionInitializer init, int p): leafParent(init.leafParent[p]), nodes(init.nodes[p].size()), color(init.leafParent[p].size(), -1){
			for (int i = 0; i < nodes.size(); i++){
				nodes[i].up = init.nodes[p][i].up;
				nodes[i].small = init.nodes[p][i].small;
				nodes[i].large = init.nodes[p][i].large;
				nodes[i].weight = init.nodes[p][i].weight;
			}
		}
		/*
		void reset(){
			totalScore = 0;
			version++;
		}
		
		void addTotal(int i){
			for (int w: leafParent[i]){
				totalZ[w].z++;
				w = totalZ[w].up;
				while (w != -1){
					Node& u = totalZ[totalZ[w].small];
					Node& v = totalZ[totalZ[w].large];
					totalZ[w].z = u.z + v.z;
					totalZ[w].z2a = u.z2a + v.z2a + u.z * v.z;
					totalZ[w].z2b = (u.z2b + v.z2b + u.z * v.z) * (1 - totalZ[w].weight);
					w = totalZ[w].up;
				}
			}
		}
		
		void rmvTotal(int i){
			for (int w: leafParent[i]){
				totalZ[w].z--;
				w = totalZ[w].up;
				while (w != -1){
					Node& u = totalZ[totalZ[w].small];
					Node& v = totalZ[totalZ[w].large];
					totalZ[w].z = u.z + v.z;
					totalZ[w].z2a = u.z2a + v.z2a + u.z * v.z;
					totalZ[w].z2b = (u.z2b + v.z2b + u.z * v.z) * (1 - totalZ[w].weight);
					w = totalZ[w].up;
				}
			}
		}
		
		void add(int x, int i){
			for (int u: leafParent[i]){
				nodes[u].update(version, totalZ[u]);
				((x == 0) ? nodes[u].z : (x == 1) ? nodes[u].x : nodes[u].y)++;
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
				((x == 0) ? nodes[u].z : (x == 1) ? nodes[u].x : nodes[u].y)--;
				int w = nodes[u].up;
				while (w != -1){
					nodes[w].update(version, totalZ[w]);
					totalScore += normal(nodes[w]);
					u = w;
					w = nodes[u].up;
				}
			}
		}
		*/
		void update(int x, int i){
			int y = color[i];
			if (x == y) return;
			for (int u: leafParent[i]){
				if (y != -1) ((y == 0) ? nodes[u].z : (y == 1) ? nodes[u].x : nodes[u].y)--;
				if (x != -1) ((x == 0) ? nodes[u].z : (x == 1) ? nodes[u].x : nodes[u].y)++;
				int w = nodes[u].up;
				while (w != -1){
					totalScore += ((nodes[w].up == -1) ? 4 : 2) * normal(nodes[w]);
					u = w;
					w = nodes[u].up;
				}
			}
			color[i] = x;
		}
		
		score_t score(){
			return totalScore;
		}
	};
	
	vector<Partition> parts;
	
	Tripartition(const TripartitionInitializer &init){
		for (int p = 0; p < init.nodes.size(); p++) parts.emplace_back(init, p);
	}
	
	void update(int x, int i){
		vector<thread> thrds;
		for (int p = 1; p < parts.size(); p++) thrds.emplace_back(&Partition::update, &parts[p], x, i);
		parts[0].update(x, i);
		for (thread &t: thrds) t.join();
	}
	
	score_t score(){
		score_t res = 0;
		for (int p = 0; p < parts.size(); p++) res += parts[p].score();
		return res;
	}
};

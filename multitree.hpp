#define OBJECTIVE_VERSION "1"

#include <vector>
#include <array>
#include <thread>
#include "threadpool.hpp"

#define SUPPORT

using namespace std;

score_t NUM_EQ_CLASSES = 0;

struct TripartitionInitializer{
	struct Node{
		int up = -1, small = -1, large = -1;
		bool dup = false;
	};
	
	vector<vector<Node> > nodes;
	vector<vector<vector<int> > > leafParent;
};

struct Tripartition{
	struct Partition{
		struct Node{
			score_t x = 0, y = 0, z = 0, tx = 0, ty = 0, tz = 0, q = 0;
			int up = -1, small = -1, large = -1; // -1 for dummy!
			bool dup = false;
		};
		
		score_t normal(Node& w){
			Node& u = nodes[w.small];
			Node& v = nodes[w.large];
			w.tx = u.tx + v.tx + u.y * v.z * (u.y + v.z - 2) + u.z * v.y * (u.z + v.y - 2);
			w.ty = u.ty + v.ty + u.x * v.z * (u.x + v.z - 2) + u.z * v.x * (u.z + v.x - 2);
			w.tz = u.tz + v.tz + u.x * v.y * (u.x + v.y - 2) + u.y * v.x * (u.y + v.x - 2);
			score_t oldQ = w.q;
			w.q = u.x * v.tx + u.y * v.ty + u.z * v.tz + v.x * u.tx + v.y * u.ty + v.z * u.tz
				+ (u.x - 1) * u.x * v.y * v.z //+ (v.x - 1) * v.x * u.y * u.z
				+ (u.y - 1) * u.y * v.x * v.z //+ (v.y - 1) * v.y * u.x * u.z
				+ (u.z - 1) * u.z * v.x * v.y //+ (v.z - 1) * v.z * u.x * u.y
				;
			return w.q - oldQ;
		}
		
		void special(Node& w){
			Node& u = nodes[w.small];
			Node& v = nodes[w.large];
			w.tx = u.tx + v.tx;
			w.ty = u.ty + v.ty;
			w.tz = u.tz + v.tz;
		}
		
		vector<vector<int> > leafParent;
		score_t totalScore = 0;
		vector<Node> nodes;
		vector<int> color;
		
		Partition(const TripartitionInitializer &init, int p): leafParent(init.leafParent[p]), nodes(init.nodes[p].size()), color(init.leafParent[p].size(), -1){
			for (int i = 0; i < nodes.size(); i++){
				nodes[i].up = init.nodes[p][i].up;
				nodes[i].small = init.nodes[p][i].small;
				nodes[i].large = init.nodes[p][i].large;
				nodes[i].dup = init.nodes[p][i].dup;
			}
		}

		void update(int x, int i){
			int y = color[i];
			if (x == y) return;
			for (int u: leafParent[i]){
				if (y != -1) ((y == 0) ? nodes[u].z : (y == 1) ? nodes[u].x : nodes[u].y)--;
				if (x != -1) ((x == 0) ? nodes[u].z : (x == 1) ? nodes[u].x : nodes[u].y)++;
				int w = nodes[u].up;
				while (w != -1 && (nodes[w].dup == false || nodes[w].large == u)){
					if (y != -1) ((y == 0) ? nodes[w].z : (y == 1) ? nodes[w].x : nodes[w].y)--;
					if (x != -1) ((x == 0) ? nodes[w].z : (x == 1) ? nodes[w].x : nodes[w].y)++;
					if (nodes[w].dup) special(nodes[w]); else totalScore += normal(nodes[w]);
					u = w;
					w = nodes[u].up;
				}
				if (w != -1){
					special(nodes[w]);
				}
			}
			color[i] = x;
		}
		
		score_t score(){
			return totalScore / 2;
		}
	};
	
	vector<Partition> parts;
	
	struct Node{
		score_t x = 0, t = 0, q = 0;
		int version = 0, up = -1, small = -1, large = -1; // -1 for dummy!
		bool dup = false;
	};
	
	Tripartition(const TripartitionInitializer &init){
		score_t s = 0;
		
		for (int p = 0; p < init.nodes.size(); p++){
			const vector<TripartitionInitializer::Node> &nodes = init.nodes[p];
			const vector<vector<int> > &leafParent = init.leafParent[p];
			vector<score_t> x(nodes.size()), t(nodes.size()), q(nodes.size());
			for (int i = 0; i < leafParent.size(); i++){
				for (int u: leafParent[i]){
					x[u]++;
					int w = nodes[u].up;
					while (w != -1 && (nodes[w].dup == false || nodes[w].large == u)){
						x[w]++;
						if (nodes[w].dup) t[w] = t[nodes[w].small] + t[nodes[w].large];
						else {
							int v = (nodes[w].small == u) ? nodes[w].large : nodes[w].small;
							t[w] = t[u] + t[v] + x[u] * x[v] * (x[u] + x[v] - 2) / 2;
							s -= q[w];
							q[w] = x[u] * t[v] + x[v] * t[u] + (x[u] - 1) * x[u] * (x[v] - 1) * x[v] / 4;
							s += q[w];
						}
						u = w;
						w = nodes[u].up;
					}
					if (w != -1) t[w] = t[nodes[w].small] + t[nodes[w].large];
				}
			}
		}
		
		NUM_EQ_CLASSES = s;
		
		for (int p = 0; p < init.nodes.size(); p++) parts.emplace_back(init, p);
	}

	void updatePart(int part, int x, int i){
		parts[part].update(x, i);
	}

	void update(int x, int i){
		vector<thread> thrds;
		for (int p = 1; p < parts.size(); p++) thrds.emplace_back(&Partition::update, &parts[p], x, i);
		parts[0].update(x, i);
		for (thread &t: thrds) t.join();
	}
	
	score_t scorePart(int part){
		return parts[part].score();
	}

	score_t score(){
		score_t res = 0;
		for (int p = 0; p < parts.size(); p++) res += parts[p].score();
		return res;
	}
};

struct Quadrupartition{
	struct Partition{
		struct Node{
			score_t a = 0, b = 0, c = 0, d = 0;
			score_t ab_c = 0, ab_d = 0, a_cd = 0, b_cd = 0;
			score_t ac_b = 0, a_bd = 0, ac_d = 0, c_bd = 0;
			score_t a_bc = 0, ad_b = 0, ad_c = 0, d_bc = 0;
			score_t ab_cd = 0, ac_bd = 0, ad_bc = 0;
			
			int up = -1, small = -1, large = -1; // -1 for dummy!
			bool dup = false;
		};
		
		array<score_t, 3> normal(Node& w){
			Node& u = nodes[w.small];
			Node& v = nodes[w.large];
			
			w.ab_c = u.ab_c + v.ab_c + u.a * u.b * v.c + u.c * v.a * v.b;
			w.ab_d = u.ab_d + v.ab_d + u.a * u.b * v.d + u.d * v.a * v.b;
			w.a_cd = u.a_cd + v.a_cd + u.c * u.d * v.a + u.a * v.c * v.d;
			w.b_cd = u.b_cd + v.b_cd + u.c * u.d * v.b + u.b * v.c * v.d;
			
			w.ac_b = u.ac_b + v.ac_b + u.a * u.c * v.b + u.b * v.a * v.c;
			w.a_bd = u.a_bd + v.a_bd + u.b * u.d * v.a + u.a * v.b * v.d;
			w.ac_d = u.ac_d + v.ac_d + u.a * u.c * v.d + u.d * v.a * v.c;
			w.c_bd = u.c_bd + v.c_bd + u.b * u.d * v.c + u.c * v.b * v.d;
			
			w.a_bc = u.a_bc + v.a_bc + u.b * u.c * v.a + u.a * v.b * v.c;
			w.ad_b = u.ad_b + v.ad_b + u.a * u.d * v.b + u.b * v.a * v.d;
			w.ad_c = u.ad_c + v.ad_c + u.a * u.d * v.c + u.c * v.a * v.d;
			w.d_bc = u.d_bc + v.d_bc + u.b * u.c * v.d + u.d * v.b * v.c;
			
			score_t ab_cd = w.ab_cd, ac_bd = w.ac_bd, ad_bc = w.ad_bc;
			w.ab_cd = u.a * u.b * v.c * v.d + u.c * u.d * v.a * v.b
					+ u.ab_c * v.d + v.ab_c * u.d
					+ u.ab_d * v.c + v.ab_d * u.c
					+ u.a_cd * v.b + v.a_cd * u.b
					+ u.b_cd * v.a + v.b_cd * u.a;
			w.ac_bd = u.a * u.c * v.b * v.d + u.b * u.d * v.a * v.c
					+ u.ac_b * v.d + v.ac_b * u.d
					+ u.a_bd * v.c + v.a_bd * u.c
					+ u.ac_d * v.b + v.ac_d * u.b
					+ u.c_bd * v.a + v.c_bd * u.a;
			w.ad_bc = u.a * u.d * v.b * v.c + u.b * u.c * v.a * v.d
					+ u.a_bc * v.d + v.a_bc * u.d
					+ u.ad_b * v.c + v.ad_b * u.c
					+ u.ad_c * v.b + v.ad_c * u.b
					+ u.d_bc * v.a + v.d_bc * u.a;
			return {w.ab_cd - ab_cd, w.ac_bd - ac_bd, w.ad_bc - ad_bc};
		}
		
		void special(Node& w){
			Node& u = nodes[w.small];
			Node& v = nodes[w.large];
			w.ab_c = u.ab_c + v.ab_c;
			w.ab_d = u.ab_d + v.ab_d;
			w.a_cd = u.a_cd + v.a_cd;
			w.b_cd = u.b_cd + v.b_cd;
			
			w.ac_b = u.ac_b + v.ac_b;
			w.a_bd = u.a_bd + v.a_bd;
			w.ac_d = u.ac_d + v.ac_d;
			w.c_bd = u.c_bd + v.c_bd;
			
			w.a_bc = u.a_bc + v.a_bc;
			w.ad_b = u.ad_b + v.ad_b;
			w.ad_c = u.ad_c + v.ad_c;
			w.d_bc = u.d_bc + v.d_bc;
		}
		
		vector<vector<int> > leafParent;
		vector<score_t> score1, score2, score3;
		score_t totalScore1 = 0, totalScore2 = 0, totalScore3 = 0, cnt[4] = {};
		vector<Node> nodes;
		vector<int> color;
		
		Partition(const TripartitionInitializer &init, int p): leafParent(init.leafParent[p]), nodes(init.nodes[p].size()), color(init.leafParent[p].size(), -1){
			for (int i = 0; i < nodes.size(); i++){
				nodes[i].up = init.nodes[p][i].up;
				nodes[i].small = init.nodes[p][i].small;
				nodes[i].large = init.nodes[p][i].large;
				nodes[i].dup = init.nodes[p][i].dup;
			}
		}
		
		void update(int x, int i){
			int y = color[i];
			if (x == y) return;
			if (y != -1) cnt[y]--;
			if (x != -1) cnt[x]++;
			for (int u: leafParent[i]){
				if (y != -1) ((y == 0) ? nodes[u].a : (y == 1) ? nodes[u].b : (y == 2) ? nodes[u].c : nodes[u].d)--;
				if (x != -1) ((x == 0) ? nodes[u].a : (x == 1) ? nodes[u].b : (x == 2) ? nodes[u].c : nodes[u].d)++;
				int w = nodes[u].up;
				while (w != -1 && (nodes[w].dup == false || nodes[w].large == u)){
					if (y != -1) ((y == 0) ? nodes[w].a : (y == 1) ? nodes[w].b : (y == 2) ? nodes[w].c : nodes[w].d)--;
					if (x != -1) ((x == 0) ? nodes[w].a : (x == 1) ? nodes[w].b : (x == 2) ? nodes[w].c : nodes[w].d)++;
					if (nodes[w].dup) special(nodes[w]);
					else {
						array<score_t, 3> s = normal(nodes[w]);
						totalScore1 += s[0]; totalScore2 += s[1]; totalScore3 += s[2];
					}
					u = w;
					w = nodes[u].up;
				}
				if (w != -1){
					special(nodes[w]);
				}
			}
			color[i] = x;
		}
		
		array<double, 3> score(){
			//cerr << cnt[0] << " " << cnt[1] << " " << cnt[2] << " " << cnt[3] << "; " << totalScore1 << " " << totalScore2 << " " << totalScore3 << "\n";
			double prod = cnt[0] * cnt[1] * cnt[2] * cnt[3];
			return {totalScore1 / prod, totalScore2 / prod, totalScore3 / prod};
		}
	};
	
	vector<Partition> parts;
	
	Quadrupartition(const TripartitionInitializer &init){
		for (int p = 0; p < init.nodes.size(); p++) parts.emplace_back(init, p);
	}
	
	void update(int x, int i){
		vector<thread> thrds;
		for (int p = 1; p < parts.size(); p++) thrds.emplace_back(&Partition::update, &parts[p], x, i);
		parts[0].update(x, i);
		for (thread &t: thrds) t.join();
	}
	
	array<double, 3> score(){
		array<double, 3> res;
		res[0] = 0;
		res[1] = 0;
		res[2] = 0;
		for (int p = 0; p < parts.size(); p++){
			auto t = parts[p].score();
			res[0] += t[0];
			res[1] += t[1];
			res[2] += t[2];
		}
		return res;
	}
};

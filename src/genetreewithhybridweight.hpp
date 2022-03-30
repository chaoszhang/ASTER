#define OBJECTIVE_VERSION "2"

#include<vector>
#include<array>
#include<thread>
#include<cmath>
#include "threadpool.hpp"

#define SUPPORT

using namespace std;

struct TripartitionInitializer{
	struct Node{
		int up = -1, small = -1, large = -1;
		score_t weight = 1, length = 1;
	};
	
	vector<int> nameCnts;
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
			score_t weight = 1, length = 1;
		};
		
		score_t normal(Node& w){
			Node& u = nodes[w.small];
			Node& v = nodes[w.large];
			//u.update(version, totalZ[w.small]);
			//v.update(version, totalZ[w.large]);
			
			w.x = (u.x + v.x) * w.length;
			w.y = (u.y + v.y) * w.length;
			w.z = (u.z + v.z) * w.length;
			w.x2a = u.x2a + v.x2a + u.x * v.x;
			w.y2a = u.y2a + v.y2a + u.y * v.y;
			w.z2a = u.z2a + v.z2a + u.z * v.z;
			w.xya = u.xya + v.xya + u.x * v.y + u.y * v.x;
			w.xza = u.xza + v.xza + u.x * v.z + u.z * v.x;
			w.yza = u.yza + v.yza + u.y * v.z + u.z * v.y;
			w.x2b = (u.x2b + v.x2b + u.x * v.x) * (1 - w.weight);
			w.y2b = (u.y2b + v.y2b + u.y * v.y) * (1 - w.weight);
			w.z2b = (u.z2b + v.z2b + u.z * v.z) * (1 - w.weight);
			w.xyb = (u.xyb + v.xyb + u.x * v.y + u.y * v.x) * (1 - w.weight);
			w.xzb = (u.xzb + v.xzb + u.x * v.z + u.z * v.x) * (1 - w.weight);
			w.yzb = (u.yzb + v.yzb + u.y * v.z + u.z * v.y) * (1 - w.weight);
			
			
			w.tx = (u.tx + v.tx + u.y * (v.z2a - v.z2b) + (u.z2a - u.z2b) * v.y + u.z * (v.y2a - v.y2b) + (u.y2a - u.y2b) * v.z) * w.length;
			w.ty = (u.ty + v.ty + u.x * (v.z2a - v.z2b) + (u.z2a - u.z2b) * v.x + u.z * (v.x2a - v.x2b) + (u.x2a - u.x2b) * v.z) * w.length;
			w.tz = (u.tz + v.tz + u.x * (v.y2a - v.y2b) + (u.y2a - u.y2b) * v.x + u.y * (v.x2a - v.x2b) + (u.x2a - u.x2b) * v.y) * w.length;
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
		
		Partition(const TripartitionInitializer &init, int p): leafParent(init.leafParent[p]), nodes(init.nodes[p].size()), color(init.leafParent[p].size(), -1){
			for (int i = 0; i < nodes.size(); i++){
				nodes[i].up = init.nodes[p][i].up;
				nodes[i].small = init.nodes[p][i].small;
				nodes[i].large = init.nodes[p][i].large;
				nodes[i].weight = init.nodes[p][i].weight;
				nodes[i].length = init.nodes[p][i].length;
			}
		}
		/*
		void reset(){
			totalScore = 0;
			version++;
		}
		
		void addTotal(int i){
			for (int w: leafParent[i]){
				totalZ[w].z += totalZ[w].length;
				w = totalZ[w].up;
				while (w != -1){
					Node& u = totalZ[totalZ[w].small];
					Node& v = totalZ[totalZ[w].large];
					totalZ[w].z = (u.z + v.z) * totalZ[w].length;
					totalZ[w].z2a = u.z2a + v.z2a + u.z * v.z;
					totalZ[w].z2b = (u.z2b + v.z2b + u.z * v.z) * (1 - totalZ[w].weight);
					w = totalZ[w].up;
				}
			}
		}
		
		void rmvTotal(int i){
			for (int w: leafParent[i]){
				totalZ[w].z -= totalZ[w].length;
				w = totalZ[w].up;
				while (w != -1){
					Node& u = totalZ[totalZ[w].small];
					Node& v = totalZ[totalZ[w].large];
					totalZ[w].z = (u.z + v.z) * totalZ[w].length;
					totalZ[w].z2a = u.z2a + v.z2a + u.z * v.z;
					totalZ[w].z2b = (u.z2b + v.z2b + u.z * v.z) * (1 - totalZ[w].weight);
					w = totalZ[w].up;
				}
			}
		}
		
		void add(int x, int i){
			for (int u: leafParent[i]){
				nodes[u].update(version, totalZ[u]);
				((x == 0) ? nodes[u].z : (x == 1) ? nodes[u].x : nodes[u].y) += nodes[u].length;
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
				((x == 0) ? nodes[u].z : (x == 1) ? nodes[u].x : nodes[u].y) -= nodes[u].length;
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
				if (y != -1) ((y == 0) ? nodes[u].z : (y == 1) ? nodes[u].x : nodes[u].y) -= nodes[u].length;
				if (x != -1) ((x == 0) ? nodes[u].z : (x == 1) ? nodes[u].x : nodes[u].y) += nodes[u].length;
				int w = nodes[u].up;
				while (w != -1){
					totalScore += normal(nodes[w]);
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
			score_t abx = 0, acx = 0, adx = 0, bcx = 0, bdx = 0, cdx = 0;
			score_t aby = 0, acy = 0, ady = 0, bcy = 0, bdy = 0, cdy = 0;
			//score_t abx = 0, cdx = 0, aby = 0, cdy = 0;
			
			score_t ab_cd = 0, ac_bd = 0, ad_bc = 0;
			
			int up = -1, small = -1, large = -1; // -1 for dummy!
			score_t weight = 1, length = 1;
		};
		
		array<score_t, 3> normal(Node& w){
			Node& u = nodes[w.small];
			Node& v = nodes[w.large];
			
			w.a = (u.a + v.a) * w.length;
			w.b = (u.b + v.b) * w.length;
			w.c = (u.c + v.c) * w.length;
			w.d = (u.d + v.d) * w.length;
			w.abx = u.abx + v.abx + u.a * v.b + v.a * u.b;
			w.acx = u.acx + v.acx + u.a * v.c + v.a * u.c;
			w.adx = u.adx + v.adx + u.a * v.d + v.a * u.d;
			w.bcx = u.bcx + v.bcx + u.b * v.c + v.b * u.c;
			w.bdx = u.bdx + v.bdx + u.b * v.d + v.b * u.d;
			w.cdx = u.cdx + v.cdx + u.c * v.d + v.c * u.d;
			w.aby = (u.aby + v.aby + u.a * v.b + v.a * u.b) * (1 - w.weight);
			w.acy = (u.acy + v.acy + u.a * v.c + v.a * u.c) * (1 - w.weight);
			w.ady = (u.ady + v.ady + u.a * v.d + v.a * u.d) * (1 - w.weight);
			w.bcy = (u.bcy + v.bcy + u.b * v.c + v.b * u.c) * (1 - w.weight);
			w.bdy = (u.bdy + v.bdy + u.b * v.d + v.b * u.d) * (1 - w.weight);
			w.cdy = (u.cdy + v.cdy + u.c * v.d + v.c * u.d) * (1 - w.weight);
			
			w.ab_c = (u.ab_c + v.ab_c + u.c * (v.abx - v.aby) + (u.abx - u.aby) * v.c) * w.length;
			w.ab_d = (u.ab_d + v.ab_d + u.d * (v.abx - v.aby) + (u.abx - u.aby) * v.d) * w.length;
			w.a_cd = (u.a_cd + v.a_cd + u.a * (v.cdx - v.cdy) + (u.cdx - u.cdy) * v.a) * w.length;
			w.b_cd = (u.b_cd + v.b_cd + u.b * (v.cdx - v.cdy) + (u.cdx - u.cdy) * v.b) * w.length;
			
			w.ac_b = (u.ac_b + v.ac_b + u.b * (v.acx - v.acy) + (u.acx - u.acy) * v.b) * w.length;
			w.ac_d = (u.ac_d + v.ac_d + u.d * (v.acx - v.acy) + (u.acx - u.acy) * v.d) * w.length;
			w.a_bd = (u.a_bd + v.a_bd + u.a * (v.bdx - v.bdy) + (u.bdx - u.bdy) * v.a) * w.length;
			w.c_bd = (u.c_bd + v.c_bd + u.c * (v.bdx - v.bdy) + (u.bdx - u.bdy) * v.c) * w.length;
			
			w.ad_b = (u.ad_b + v.ad_b + u.b * (v.adx - v.ady) + (u.adx - u.ady) * v.b) * w.length;
			w.ad_c = (u.ad_c + v.ad_c + u.c * (v.adx - v.ady) + (u.adx - u.ady) * v.c) * w.length;
			w.a_bc = (u.a_bc + v.a_bc + u.a * (v.bcx - v.bcy) + (u.bcx - u.bcy) * v.a) * w.length;
			w.d_bc = (u.d_bc + v.d_bc + u.d * (v.bcx - v.bcy) + (u.bcx - u.bcy) * v.d) * w.length;
			
			score_t ab_cd = w.ab_cd, ac_bd = w.ac_bd, ad_bc = w.ad_bc;
			w.ab_cd = u.abx * v.cdx - u.aby * v.cdy + v.abx * u.cdx - v.aby * u.cdy
					+ u.ab_c * v.d + v.ab_c * u.d
					+ u.ab_d * v.c + v.ab_d * u.c
					+ u.a_cd * v.b + v.a_cd * u.b
					+ u.b_cd * v.a + v.b_cd * u.a;
			w.ac_bd = u.acx * v.bdx - u.acy * v.bdy + v.acx * u.bdx - v.acy * u.bdy
					+ u.ac_b * v.d + v.ac_b * u.d
					+ u.a_bd * v.c + v.a_bd * u.c
					+ u.ac_d * v.b + v.ac_d * u.b
					+ u.c_bd * v.a + v.c_bd * u.a;
			w.ad_bc = u.adx * v.bcx - u.ady * v.bcy + v.adx * u.bcx - v.ady * u.bcy
					+ u.a_bc * v.d + v.a_bc * u.d
					+ u.ad_b * v.c + v.ad_b * u.c
					+ u.ad_c * v.b + v.ad_c * u.b
					+ u.d_bc * v.a + v.d_bc * u.a;
			return {w.ab_cd - ab_cd, w.ac_bd - ac_bd, w.ad_bc - ad_bc};
		}
		
		vector<vector<int> > leafParent;
		score_t totalScore1 = 0, totalScore2 = 0, totalScore3 = 0, cnt[4] = {}, totalSqrScore = 0, gtCnt = 0;
		vector<Node> nodes;
		vector<score_t> rootScores;
		vector<int> color;
		vector<int> nameCnts;
		
		Partition(TripartitionInitializer init, int p): leafParent(init.leafParent[p]), nodes(init.nodes[p].size()), rootScores(init.nodes[p].size()), color(init.leafParent[p].size(), -1), nameCnts(init.nameCnts){
			for (int i = 0; i < nodes.size(); i++){
				nodes[i].up = init.nodes[p][i].up;
				nodes[i].small = init.nodes[p][i].small;
				nodes[i].large = init.nodes[p][i].large;
				nodes[i].weight = init.nodes[p][i].weight;
				nodes[i].length = init.nodes[p][i].length;
			}
		}
		
		void update(int x, int i){
			int y = color[i];
			if (x == y) return;
			if (y != -1) cnt[y] -= nameCnts[i];
			if (x != -1) cnt[x] += nameCnts[i];
			for (int u: leafParent[i]){
				score_t s1 = 0, s2 = 0, s3 = 0;
				if (y != -1) ((y == 0) ? nodes[u].a : (y == 1) ? nodes[u].b : (y == 2) ? nodes[u].c : nodes[u].d) -= nodes[u].length;
				if (x != -1) ((x == 0) ? nodes[u].a : (x == 1) ? nodes[u].b : (x == 2) ? nodes[u].c : nodes[u].d) += nodes[u].length;
				int w = nodes[u].up;
				while (w != -1){
					array<score_t, 3> s = normal(nodes[w]);
					s1 += s[0]; s2 += s[1]; s3 += s[2];
					u = w;
					w = nodes[u].up;
				}
				totalScore1 += s1;
				totalScore2 += s2;
				totalScore3 += s3;
				totalSqrScore -= rootScores[u] * rootScores[u];
				if (rootScores[u] > 1e-9) gtCnt--;
				rootScores[u] += s1 + s2 + s3;
				totalSqrScore += rootScores[u] * rootScores[u];
				if (rootScores[u] > 1e-9) gtCnt++;
			}
			color[i] = x;
		}
		
		array<score_t, 5> score(){
			return {totalScore1, totalScore2, totalScore3, totalSqrScore, gtCnt};
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
		double sqrsum = 0;
		double cnt = 0;
		for (int p = 0; p < parts.size(); p++){
			auto t = parts[p].score();
			res[0] += t[0];
			res[1] += t[1];
			res[2] += t[2];
			sqrsum += t[3];
			cnt += t[4];
		}
		res[0] = (cnt < 0.5) ? 0 : res[0] * sqrt(cnt / sqrsum);
		res[1] = (cnt < 0.5) ? 0 : res[1] * sqrt(cnt / sqrsum);
		res[2] = (cnt < 0.5) ? 0 : res[2] * sqrt(cnt / sqrsum);
		return res;
	}
};

#define OBJECTIVE_VERSION "4"

/* CHANGE LOG
 * 4. adding support for CASTLES
 * 3. code refactoring
 * 2. add customized annotation for branch lengths
 * 1. replace genetree.hpp as the default objective for astral.cpp
 */

#include<vector>
#include<array>
#include<thread>
#include "threadpool.hpp"

#define SUPPORT
#define CUSTOMIZED_ANNOTATION

using namespace std;

struct CustomizedAnnotation{
	struct Quadrupartition{
		score_t quartetCnt = 0;
		length_t sumInternalLength = 0;
		length_t sumLengthA = 0;
		length_t sumLengthB = 0;
		length_t sumLengthC = 0;
		length_t sumLengthD = 0;

		void numericalBalancing(){
			if (quartetCnt == 0){
				sumInternalLength = 0;
				sumLengthA = 0;
				sumLengthB = 0;
				sumLengthC = 0;
				sumLengthD = 0;
			}
		}
		
		Quadrupartition operator+ (const Quadrupartition& o) const{
			Quadrupartition r;
			r.quartetCnt = quartetCnt + o.quartetCnt;
			r.sumInternalLength = sumInternalLength + o.sumInternalLength;
			r.sumLengthA = sumLengthA + o.sumLengthA;
			r.sumLengthB = sumLengthB + o.sumLengthB;
			r.sumLengthC = sumLengthC + o.sumLengthC;
			r.sumLengthD = sumLengthD + o.sumLengthD;
			return r;
		}

		Quadrupartition& operator+= (const Quadrupartition& o){
			quartetCnt += o.quartetCnt;
			sumInternalLength += o.sumInternalLength;
			sumLengthA += o.sumLengthA;
			sumLengthB += o.sumLengthB;
			sumLengthC += o.sumLengthC;
			sumLengthD += o.sumLengthD;
			return *this;
		}

		Quadrupartition operator- (const Quadrupartition& o) const{
			Quadrupartition r;
			r.quartetCnt = quartetCnt - o.quartetCnt;
			r.sumInternalLength = sumInternalLength - o.sumInternalLength;
			r.sumLengthA = sumLengthA - o.sumLengthA;
			r.sumLengthB = sumLengthB - o.sumLengthB;
			r.sumLengthC = sumLengthC - o.sumLengthC;
			r.sumLengthD = sumLengthD - o.sumLengthD;
			return r;
		}
	} ab_cd, ac_bd, ad_bc;
	// a: left child, b: right child, c: sibling, d: outgroup

	CustomizedAnnotation operator+ (const CustomizedAnnotation& o) const{
		CustomizedAnnotation r;
		r.ab_cd = ab_cd + o.ab_cd;
		r.ac_bd = ac_bd + o.ac_bd;
		r.ad_bc = ad_bc + o.ad_bc;
		return r;
	}

	CustomizedAnnotation& operator+= (const CustomizedAnnotation& o){
		ab_cd += o.ab_cd;
		ac_bd += o.ac_bd;
		ad_bc += o.ad_bc;
		return *this;
	}

	CustomizedAnnotation operator- (const CustomizedAnnotation& o) const{
		CustomizedAnnotation r;
		r.ab_cd = ab_cd - o.ab_cd;
		r.ac_bd = ac_bd - o.ac_bd;
		r.ad_bc = ad_bc - o.ad_bc;
		return r;
	}
};

struct TripartitionInitializer{
	struct Node{
		int up = -1, small = -1, large = -1;
		bool isGhostBranch = 1;
		length_t length = 0;
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
			bool isGhostBranch = 1;
		};
		
		score_t normal(Node& w){
			Node& u = nodes[w.small];
			Node& v = nodes[w.large];
			
			w.x = u.x + v.x;
			w.y = u.y + v.y;
			w.z = u.z + v.z;
			w.x2a = u.x2a + v.x2a + u.x * v.x;
			w.y2a = u.y2a + v.y2a + u.y * v.y;
			w.z2a = u.z2a + v.z2a + u.z * v.z;
			w.xya = u.xya + v.xya + u.x * v.y + u.y * v.x;
			w.xza = u.xza + v.xza + u.x * v.z + u.z * v.x;
			w.yza = u.yza + v.yza + u.y * v.z + u.z * v.y;
			if (w.isGhostBranch) {
				w.x2b = u.x2b + v.x2b + u.x * v.x;
				w.y2b = u.y2b + v.y2b + u.y * v.y;
				w.z2b = u.z2b + v.z2b + u.z * v.z;
				w.xyb = u.xyb + v.xyb + u.x * v.y + u.y * v.x;
				w.xzb = u.xzb + v.xzb + u.x * v.z + u.z * v.x;
				w.yzb = u.yzb + v.yzb + u.y * v.z + u.z * v.y;
			}
			else w.x2b = w.y2b = w.z2b = w.xyb = w.xzb = w.yzb = 0;
			w.tx = u.tx + v.tx + u.y * (v.z2a - v.z2b) + (u.z2a - u.z2b) * v.y + u.z * (v.y2a - v.y2b) + (u.y2a - u.y2b) * v.z;
			w.ty = u.ty + v.ty + u.x * (v.z2a - v.z2b) + (u.z2a - u.z2b) * v.x + u.z * (v.x2a - v.x2b) + (u.x2a - u.x2b) * v.z;
			w.tz = u.tz + v.tz + u.x * (v.y2a - v.y2b) + (u.y2a - u.y2b) * v.x + u.y * (v.x2a - v.x2b) + (u.x2a - u.x2b) * v.y;
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
				nodes[i].isGhostBranch = init.nodes[p][i].isGhostBranch;
			}
		}
		
		void update(int x, int i){
			int y = color[i];
			if (x == y) return;
			for (int u: leafParent[i]){
				if (y != -1) ((y == 0) ? nodes[u].z : (y == 1) ? nodes[u].x : nodes[u].y)--;
				if (x != -1) ((x == 0) ? nodes[u].z : (x == 1) ? nodes[u].x : nodes[u].y)++;
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
			struct Topology{
				score_t pq_r = 0, pq_s = 0, p_rs = 0, q_rs = 0;
				score_t pqx = 0, rsx = 0, pqy = 0, rsy = 0;
				score_t pq_rs = 0;
				length_t Pqx = 0, Qpx = 0, Rsx = 0, Srx = 0;
				length_t Pqy = 0, Qpy = 0, Rsy = 0, Sry = 0;
				length_t pqX = 0, rsX = 0;
				length_t pqY = 0, rsY = 0; 
				length_t Pq_r = 0, Pq_s = 0, P_rs = 0;
				length_t Qp_r = 0, Qp_s = 0, Q_rs = 0;
				length_t Rs_q = 0, Rs_p = 0, R_pq = 0;
				length_t Sr_q = 0, Sr_p = 0, S_pq = 0;
				length_t pqZr = 0, pqZs = 0, pZrs = 0, qZrs = 0;
				length_t Pq_rs = 0, Qp_rs = 0, Rs_pq = 0, Sr_pq = 0;
			} ab_cd, ac_bd, ad_bc;

			score_t a = 0, b = 0, c = 0, d = 0;
			length_t A = 0, B = 0, C = 0, D = 0;
			int up = -1, small = -1, large = -1; // -1 for dummy!
			bool isGhostBranch = 1;
			length_t length = 0;
			CustomizedAnnotation annot;
		};

		template<int T> void extended(Node& w, Node& u, Node& v){
			Node::Topology& U = ((T == 0) ? u.ab_cd : (T == 1) ? u.ac_bd : u.ad_bc);
			Node::Topology& V = ((T == 0) ? v.ab_cd : (T == 1) ? v.ac_bd : v.ad_bc);
			Node::Topology& W = ((T == 0) ? w.ab_cd : (T == 1) ? w.ac_bd : w.ad_bc);
			score_t up = u.a, vp = v.a;
			score_t uq = ((T == 0) ? u.b : (T == 1) ? u.c : u.d);
			score_t ur = ((T == 0) ? u.c : (T == 1) ? u.b : u.b);
			score_t us = ((T == 0) ? u.d : (T == 1) ? u.d : u.c);
			score_t vq = ((T == 0) ? v.b : (T == 1) ? v.c : v.d);
			score_t vr = ((T == 0) ? v.c : (T == 1) ? v.b : v.b);
			score_t vs = ((T == 0) ? v.d : (T == 1) ? v.d : v.c);
			length_t uP = u.A, vP = v.A;
			length_t uQ = ((T == 0) ? u.B : (T == 1) ? u.C : u.D);
			length_t uR = ((T == 0) ? u.C : (T == 1) ? u.B : u.B);
			length_t uS = ((T == 0) ? u.D : (T == 1) ? u.D : u.C);
			length_t vQ = ((T == 0) ? v.B : (T == 1) ? v.C : v.D);
			length_t vR = ((T == 0) ? v.C : (T == 1) ? v.B : v.B);
			length_t vS = ((T == 0) ? v.D : (T == 1) ? v.D : v.C);
			CustomizedAnnotation::Quadrupartition& annot = ((T == 0) ? w.annot.ab_cd : (T == 1) ? w.annot.ac_bd : w.annot.ad_bc);
			
			W.pqx = U.pqx + V.pqx + up * vq + vp * uq;
			W.rsx = U.rsx + V.rsx + ur * vs + vr * us;

			if (w.isGhostBranch) {
				W.pqy = U.pqy + V.pqy + up * vq + vp * uq;
				W.rsy = U.rsy + V.rsy + ur * vs + vr * us;
			}
			else {
				W.pqy = W.rsy = 0;
			}

			W.pq_r = U.pq_r + V.pq_r + ur * (V.pqx - V.pqy) + (U.pqx - U.pqy) * vr;
			W.pq_s = U.pq_s + V.pq_s + us * (V.pqx - V.pqy) + (U.pqx - U.pqy) * vs;
			W.p_rs = U.p_rs + V.p_rs + up * (V.rsx - V.rsy) + (U.rsx - U.rsy) * vp;
			W.q_rs = U.q_rs + V.q_rs + uq * (V.rsx - V.rsy) + (U.rsx - U.rsy) * vq;
			
			W.pq_rs = U.pqx * V.rsx - U.pqy * V.rsy + V.pqx * U.rsx - V.pqy * U.rsy
					+ U.pq_r * vs + V.pq_r * us
					+ U.pq_s * vr + V.pq_s * ur
					+ U.p_rs * vq + V.p_rs * uq
					+ U.q_rs * vp + V.q_rs * up;
			
			W.Pqx = U.Pqx + V.Pqx + uP * vq + vP * uq;
			W.Qpx = U.Qpx + V.Qpx + uQ * vp + vQ * up;
			W.Rsx = U.Rsx + V.Rsx + uR * vs + vR * us;
			W.Srx = U.Srx + V.Srx + uS * vr + vS * ur;

			W.pqX = U.pqX + V.pqX + W.pqx * w.length;
			W.rsX = U.rsX + V.rsX + W.rsx * w.length;

			if (w.isGhostBranch) {
				W.Pqy = U.Pqy + V.Pqy + uP * vq + vP * uq;
				W.Qpy = U.Qpy + V.Qpy + uQ * vp + vQ * up;
				W.Rsy = U.Rsy + V.Rsy + uR * vs + vR * us;
				W.Sry = U.Sry + V.Sry + uS * vr + vS * ur;
				W.pqY = U.pqY + V.pqY + W.pqy * w.length;
				W.rsY = U.rsY + V.rsY + W.rsy * w.length;
			}
			else {
				W.Pqy = W.Qpy = W.Rsy = W.Sry = 0;
				W.pqY = W.rsY = 0;
			}

			W.Pq_r = U.Pq_r + V.Pq_r + (U.Pqx - U.Pqy) * vr + (V.Pqx - V.Pqy) * ur;
			W.Pq_s = U.Pq_s + V.Pq_s + (U.Pqx - U.Pqy) * vs + (V.Pqx - V.Pqy) * us;
			W.Qp_r = U.Qp_r + V.Qp_r + (U.Qpx - U.Qpy) * vr + (V.Qpx - V.Qpy) * ur;
			W.Qp_s = U.Qp_s + V.Qp_s + (U.Qpx - U.Qpy) * vs + (V.Qpx - V.Qpy) * us;
			W.Rs_p = U.Rs_p + V.Rs_p + (U.Rsx - U.Rsy) * vp + (V.Rsx - V.Rsy) * up;
			W.Rs_q = U.Rs_q + V.Rs_q + (U.Rsx - U.Rsy) * vq + (V.Rsx - V.Rsy) * uq;
			W.Sr_p = U.Sr_p + V.Sr_p + (U.Srx - U.Sry) * vp + (V.Srx - V.Sry) * up;
			W.Sr_q = U.Sr_q + V.Sr_q + (U.Srx - U.Sry) * vq + (V.Srx - V.Sry) * uq;
			
			W.P_rs = U.P_rs + V.P_rs + uP * (V.rsx - V.rsy) + vP * (U.rsx - U.rsy);
			W.Q_rs = U.Q_rs + V.Q_rs + uQ * (V.rsx - V.rsy) + vQ * (U.rsx - U.rsy);
			W.R_pq = U.R_pq + V.R_pq + uR * (V.pqx - V.pqy) + vR * (U.pqx - U.pqy);
			W.S_pq = U.S_pq + V.S_pq + uS * (V.pqx - V.pqy) + vS * (U.pqx - U.pqy);
			
			W.pqZr = U.pqZr + V.pqZr + (U.pqX - U.pqY) * vr + (V.pqX - V.pqY) * ur;
			W.pqZs = U.pqZs + V.pqZs + (U.pqX - U.pqY) * vs + (V.pqX - V.pqY) * us;
			W.pZrs = U.pZrs + V.pZrs + (U.rsX - U.rsY) * vp + (V.rsX - V.rsY) * up;
			W.qZrs = U.qZrs + V.qZrs + (U.rsX - U.rsY) * vq + (V.rsX - V.rsY) * uq;

			W.Pq_rs = U.Pq_rs + V.Pq_rs + W.q_rs * w.length;
			W.Qp_rs = U.Qp_rs + V.Qp_rs + W.p_rs * w.length;
			W.Rs_pq = U.Rs_pq + V.Rs_pq + W.pq_s * w.length;
			W.Sr_pq = U.Sr_pq + V.Sr_pq + W.pq_r * w.length;

			annot.quartetCnt = W.pq_rs;
			annot.sumLengthA = U.Pq_r * vs + V.Pq_r * us + U.Pq_s * vr + V.Pq_s * ur + U.P_rs * vq + V.P_rs * uq + uP * V.q_rs + vP * U.q_rs + up * V.Pq_rs + vp * U.Pq_rs + U.Pqx * V.rsx - U.Pqy * V.rsy + V.Pqx * U.rsx - V.Pqy * U.rsy;
			annot.sumLengthB = U.Qp_r * vs + V.Qp_r * us + U.Qp_s * vr + V.Qp_s * ur + U.Q_rs * vp + V.Q_rs * up + uQ * V.p_rs + vQ * U.p_rs + uq * V.Qp_rs + vq * U.Qp_rs + U.Qpx * V.rsx - U.Qpy * V.rsy + V.Qpx * U.rsx - V.Qpy * U.rsy;
			annot.sumLengthC = U.Rs_p * vq + V.Rs_p * uq + U.Rs_q * vp + V.Rs_q * up + U.R_pq * vs + V.R_pq * us + uR * V.pq_s + vR * U.pq_s + ur * V.Rs_pq + vr * U.Rs_pq + U.Rsx * V.pqx - U.Rsy * V.pqy + V.Rsx * U.pqx - V.Rsy * U.pqy;
			annot.sumLengthD = U.Sr_p * vq + V.Sr_p * uq + U.Sr_q * vp + V.Sr_q * up + U.S_pq * vr + V.S_pq * ur + uS * V.pq_r + vS * U.pq_r + us * V.Sr_pq + vs * U.Sr_pq + U.Srx * V.pqx - U.Sry * V.pqy + V.Srx * U.pqx - V.Sry * U.pqy;
			annot.sumInternalLength = U.pqZr * vs + V.pqZr * us + U.pqZs * vr + V.pqZs * ur + U.pZrs * vq + V.pZrs * uq + up * V.qZrs + vp * U.qZrs
									+ U.pqX * V.rsx - U.pqY * V.rsy + U.pqx * V.rsX - U.pqy * V.rsY + V.pqX * U.rsx - V.pqY * U.rsy + V.pqx * U.rsX - V.pqy * U.rsY;
		}

		CustomizedAnnotation extended(Node& w){
			CustomizedAnnotation old = w.annot;

			Node& u = nodes[w.small];
			Node& v = nodes[w.large];

			w.a = u.a + v.a;
			w.b = u.b + v.b;
			w.c = u.c + v.c;
			w.d = u.d + v.d;
			
			w.A = u.A + v.A + w.a * w.length;
			w.B = u.B + v.B + w.b * w.length;
			w.C = u.C + v.C + w.c * w.length;
			w.D = u.D + v.D + w.d * w.length;

			extended<0>(w, u, v);
			extended<1>(w, u, v);
			extended<2>(w, u, v);

			return w.annot - old;
		}
		
		vector<vector<int> > leafParent;
		vector<score_t> score1, score2, score3, prod;
		double totalScore1 = 0, totalScore2 = 0, totalScore3 = 0;
		CustomizedAnnotation annot;
		vector<Node> nodes;
		vector<int> color;
		
		Partition(TripartitionInitializer init, int p): leafParent(init.leafParent[p]), nodes(init.nodes[p].size()), 
				prod(init.nodes[p].size()), score1(init.nodes[p].size()), score2(init.nodes[p].size()), score3(init.nodes[p].size()), color(init.leafParent[p].size(), -1){
			for (int i = 0; i < nodes.size(); i++){
				nodes[i].up = init.nodes[p][i].up;
				nodes[i].small = init.nodes[p][i].small;
				nodes[i].large = init.nodes[p][i].large;
				nodes[i].isGhostBranch = init.nodes[p][i].isGhostBranch;
				nodes[i].length = init.nodes[p][i].length;
			}
		}
		
		void update(int x, int i){
			int y = color[i];
			if (x == y) return;
			for (int u: leafParent[i]){
				score_t s1 = 0, s2 = 0, s3 = 0;
				if (y != -1) ((y == 0) ? nodes[u].a : (y == 1) ? nodes[u].b : (y == 2) ? nodes[u].c : nodes[u].d)--;
				if (x != -1) ((x == 0) ? nodes[u].a : (x == 1) ? nodes[u].b : (x == 2) ? nodes[u].c : nodes[u].d)++;
				if (y != -1) ((y == 0) ? nodes[u].A : (y == 1) ? nodes[u].B : (y == 2) ? nodes[u].C : nodes[u].D) -= nodes[u].length;
				if (x != -1) ((x == 0) ? nodes[u].A : (x == 1) ? nodes[u].B : (x == 2) ? nodes[u].C : nodes[u].D) += nodes[u].length;
				int w = nodes[u].up;
				while (w != -1){
					CustomizedAnnotation diff = extended(nodes[w]);
					annot += diff;
					s1 += diff.ab_cd.quartetCnt;
					s2 += diff.ac_bd.quartetCnt;
					s3 += diff.ad_bc.quartetCnt;
					u = w;
					w = nodes[u].up;
				}
				score_t t = score1[u] + score2[u] + score3[u];
				length_t dt = t;
				if (t > 0){
					totalScore1 -= score1[u] / dt;
					totalScore2 -= score2[u] / dt;
					totalScore3 -= score3[u] / dt;
				}
				score1[u] += s1; score2[u] += s2; score3[u] += s3;
				t = score1[u] + score2[u] + score3[u];
				dt = t;
				if (t > 0){
					totalScore1 += score1[u] / dt;
					totalScore2 += score2[u] / dt;
					totalScore3 += score3[u] / dt;
				}
			}
			color[i] = x;
		}
		
		array<double, 3> score(){
			return {(double) totalScore1, (double) totalScore2, (double) totalScore3};
		}

		CustomizedAnnotation annotate(){
			annot.ab_cd.numericalBalancing();
			annot.ac_bd.numericalBalancing();
			annot.ad_bc.numericalBalancing();
			return annot;
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

	CustomizedAnnotation annotate(){
		CustomizedAnnotation res;
		for (int p = 0; p < parts.size(); p++){
			res += parts[p].annotate();
		}
		return res;
	}
};

#define OBJECTIVE_VERSION "3"

/* CHANGE LOG
 * 3. add CASTLES-Pro branch length calculations
 * 2. add customized annotation for branch lengths
 */

#include <vector>
#include <array>
#include <thread>
#include "threadpool.hpp"

#define SUPPORT
#define CUSTOMIZED_ANNOTATION

using namespace std;

score_t NUM_EQ_CLASSES = 0;

struct CustomizedAnnotation{
	struct Quadrupartition{
		length_t quartetCnt = 0;
		length_t sumInternalLength = 0;
		length_t sumLengthA = 0;
		length_t sumLengthB = 0;
		length_t sumLengthC = 0;
		length_t sumLengthD = 0;

		void numericalBalancing(){
			if (quartetCnt < 1e-6){
				quartetCnt = 0;
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
		bool dup = false;
		length_t length;
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
			struct Topology{
				length_t pq_r = 0, pq_s = 0, p_rs = 0, q_rs = 0;
				length_t pq = 0, rs = 0;
				length_t pq_rs = 0;
				length_t Pq = 0, Qp = 0, Rs = 0, Sr = 0;
				length_t pqX = 0, rsX = 0;
				length_t Pq_r = 0, Pq_s = 0, P_rs = 0;
				length_t Qp_r = 0, Qp_s = 0, Q_rs = 0;
				length_t Rs_q = 0, Rs_p = 0, R_pq = 0;
				length_t Sr_q = 0, Sr_p = 0, S_pq = 0;
				length_t pqXr = 0, pqXs = 0, pXrs = 0, qXrs = 0;
				length_t Pq_rs = 0, Qp_rs = 0, Rs_pq = 0, Sr_pq = 0;
			} T_ab_cd, T_ac_bd, T_ad_bc;

			length_t aa = 0, bb = 0, cc = 0, dd = 0;
			length_t A = 0, B = 0, C = 0, D = 0;
			// bool isGhostBranch = 1;
			length_t length = 0;
			CustomizedAnnotation annot;

			score_t a = 0, b = 0, c = 0, d = 0;
			score_t ab_c = 0, ab_d = 0, a_cd = 0, b_cd = 0;
			score_t ac_b = 0, a_bd = 0, ac_d = 0, c_bd = 0;
			score_t a_bc = 0, ad_b = 0, ad_c = 0, d_bc = 0;
			score_t ab_cd = 0, ac_bd = 0, ad_bc = 0;
			
			int up = -1, small = -1, large = -1; // -1 for dummy!
			bool dup = false;
		};
		
		template<int T> void extended(Node& w, Node& u, Node& v){
			Node::Topology& U = ((T == 0) ? u.T_ab_cd : (T == 1) ? u.T_ac_bd : u.T_ad_bc);
			Node::Topology& V = ((T == 0) ? v.T_ab_cd : (T == 1) ? v.T_ac_bd : v.T_ad_bc);
			Node::Topology& W = ((T == 0) ? w.T_ab_cd : (T == 1) ? w.T_ac_bd : w.T_ad_bc);
			score_t up = u.aa, vp = v.aa;
			score_t uq = ((T == 0) ? u.bb : (T == 1) ? u.cc : u.dd);
			score_t ur = ((T == 0) ? u.cc : (T == 1) ? u.bb : u.bb);
			score_t us = ((T == 0) ? u.dd : (T == 1) ? u.dd : u.cc);
			score_t vq = ((T == 0) ? v.bb : (T == 1) ? v.cc : v.dd);
			score_t vr = ((T == 0) ? v.cc : (T == 1) ? v.bb : v.bb);
			score_t vs = ((T == 0) ? v.dd : (T == 1) ? v.dd : v.cc);
			length_t uP = u.A, vP = v.A;
			length_t uQ = ((T == 0) ? u.B : (T == 1) ? u.C : u.D);
			length_t uR = ((T == 0) ? u.C : (T == 1) ? u.B : u.B);
			length_t uS = ((T == 0) ? u.D : (T == 1) ? u.D : u.C);
			length_t vQ = ((T == 0) ? v.B : (T == 1) ? v.C : v.D);
			length_t vR = ((T == 0) ? v.C : (T == 1) ? v.B : v.B);
			length_t vS = ((T == 0) ? v.D : (T == 1) ? v.D : v.C);
			CustomizedAnnotation::Quadrupartition& annot = ((T == 0) ? w.annot.ab_cd : (T == 1) ? w.annot.ac_bd : w.annot.ad_bc);

			if (w.dup == false){
				W.pq = U.pq + V.pq + up * vq + vp * uq;
				W.rs = U.rs + V.rs + ur * vs + vr * us;

				W.pq_r = U.pq_r + V.pq_r + ur * V.pq + U.pq * vr;
				W.pq_s = U.pq_s + V.pq_s + us * V.pq + U.pq * vs;
				W.p_rs = U.p_rs + V.p_rs + up * V.rs + U.rs * vp;
				W.q_rs = U.q_rs + V.q_rs + uq * V.rs + U.rs * vq;
				
				W.pq_rs = U.pq * V.rs + V.pq * U.rs
						+ U.pq_r * vs + V.pq_r * us
						+ U.pq_s * vr + V.pq_s * ur
						+ U.p_rs * vq + V.p_rs * uq
						+ U.q_rs * vp + V.q_rs * up;
				
				W.Pq = U.Pq + V.Pq + uP * vq + vP * uq;
				W.Qp = U.Qp + V.Qp + uQ * vp + vQ * up;
				W.Rs = U.Rs + V.Rs + uR * vs + vR * us;
				W.Sr = U.Sr + V.Sr + uS * vr + vS * ur;

				W.pqX = U.pqX + V.pqX + W.pq * w.length;
				W.rsX = U.rsX + V.rsX + W.rs * w.length;

				W.Pq_r = U.Pq_r + V.Pq_r + U.Pq * vr + V.Pq * ur;
				W.Pq_s = U.Pq_s + V.Pq_s + U.Pq * vs + V.Pq * us;
				W.Qp_r = U.Qp_r + V.Qp_r + U.Qp * vr + V.Qp * ur;
				W.Qp_s = U.Qp_s + V.Qp_s + U.Qp * vs + V.Qp * us;
				W.Rs_p = U.Rs_p + V.Rs_p + U.Rs * vp + V.Rs * up;
				W.Rs_q = U.Rs_q + V.Rs_q + U.Rs * vq + V.Rs * uq;
				W.Sr_p = U.Sr_p + V.Sr_p + U.Sr * vp + V.Sr * up;
				W.Sr_q = U.Sr_q + V.Sr_q + U.Sr * vq + V.Sr * uq;
				
				W.P_rs = U.P_rs + V.P_rs + uP * V.rs + vP * U.rs;
				W.Q_rs = U.Q_rs + V.Q_rs + uQ * V.rs + vQ * U.rs;
				W.R_pq = U.R_pq + V.R_pq + uR * V.pq + vR * U.pq;
				W.S_pq = U.S_pq + V.S_pq + uS * V.pq + vS * U.pq;
				
				W.pqXr = U.pqXr + V.pqXr + U.pqX * vr + V.pqX * ur;
				W.pqXs = U.pqXs + V.pqXs + U.pqX * vs + V.pqX * us;
				W.pXrs = U.pXrs + V.pXrs + U.rsX * vp + V.rsX * up;
				W.qXrs = U.qXrs + V.qXrs + U.rsX * vq + V.rsX * uq;

				W.Pq_rs = U.Pq_rs + V.Pq_rs + W.q_rs * w.length;
				W.Qp_rs = U.Qp_rs + V.Qp_rs + W.p_rs * w.length;
				W.Rs_pq = U.Rs_pq + V.Rs_pq + W.pq_s * w.length;
				W.Sr_pq = U.Sr_pq + V.Sr_pq + W.pq_r * w.length;

				annot.quartetCnt = W.pq_rs;
				annot.sumLengthA = U.Pq_r * vs + V.Pq_r * us + U.Pq_s * vr + V.Pq_s * ur + U.P_rs * vq + V.P_rs * uq + uP * V.q_rs + vP * U.q_rs + up * V.Pq_rs + vp * U.Pq_rs + U.Pq * V.rs + V.Pq * U.rs;
				annot.sumLengthB = U.Qp_r * vs + V.Qp_r * us + U.Qp_s * vr + V.Qp_s * ur + U.Q_rs * vp + V.Q_rs * up + uQ * V.p_rs + vQ * U.p_rs + uq * V.Qp_rs + vq * U.Qp_rs + U.Qp * V.rs + V.Qp * U.rs;
				annot.sumLengthC = U.Rs_p * vq + V.Rs_p * uq + U.Rs_q * vp + V.Rs_q * up + U.R_pq * vs + V.R_pq * us + uR * V.pq_s + vR * U.pq_s + ur * V.Rs_pq + vr * U.Rs_pq + U.Rs * V.pq + V.Rs * U.pq;
				annot.sumLengthD = U.Sr_p * vq + V.Sr_p * uq + U.Sr_q * vp + V.Sr_q * up + U.S_pq * vr + V.S_pq * ur + uS * V.pq_r + vS * U.pq_r + us * V.Sr_pq + vs * U.Sr_pq + U.Sr * V.pq + V.Sr * U.pq;
				annot.sumInternalLength = U.pqXr * vs + V.pqXr * us + U.pqXs * vr + V.pqXs * ur + U.pXrs * vq + V.pXrs * uq + up * V.qXrs + vp * U.qXrs
										+ U.pqX * V.rs + U.pq * V.rsX + V.pqX * U.rs + V.pq * U.rsX;
			}
			else {
				W.pq = (U.pq + V.pq) / 2;
				W.rs = (U.rs + V.rs) / 2;

				W.pq_r = (U.pq_r + V.pq_r) / 2;
				W.pq_s = (U.pq_s + V.pq_s) / 2;
				W.p_rs = (U.p_rs + V.p_rs) / 2;
				W.q_rs = (U.q_rs + V.q_rs) / 2;
				
				W.pq_rs = 0;
				
				W.Pq = (U.Pq + V.Pq) / 2;
				W.Qp = (U.Qp + V.Qp) / 2;
				W.Rs = (U.Rs + V.Rs) / 2;
				W.Sr = (U.Sr + V.Sr) / 2;

				W.pqX = (U.pqX + V.pqX) / 2 + W.pq * w.length;
				W.rsX = (U.rsX + V.rsX) / 2 + W.rs * w.length;

				W.Pq_r = (U.Pq_r + V.Pq_r) / 2;
				W.Pq_s = (U.Pq_s + V.Pq_s) / 2;
				W.Qp_r = (U.Qp_r + V.Qp_r) / 2;
				W.Qp_s = (U.Qp_s + V.Qp_s) / 2;
				W.Rs_p = (U.Rs_p + V.Rs_p) / 2;
				W.Rs_q = (U.Rs_q + V.Rs_q) / 2;
				W.Sr_p = (U.Sr_p + V.Sr_p) / 2;
				W.Sr_q = (U.Sr_q + V.Sr_q) / 2;
				
				W.P_rs = (U.P_rs + V.P_rs) / 2;
				W.Q_rs = (U.Q_rs + V.Q_rs) / 2;
				W.R_pq = (U.R_pq + V.R_pq) / 2;
				W.S_pq = (U.S_pq + V.S_pq) / 2;
				
				W.pqXr = (U.pqXr + V.pqXr) / 2;
				W.pqXs = (U.pqXs + V.pqXs) / 2;
				W.pXrs = (U.pXrs + V.pXrs) / 2;
				W.qXrs = (U.qXrs + V.qXrs) / 2;

				W.Pq_rs = (U.Pq_rs + V.Pq_rs) / 2 + W.q_rs * w.length;
				W.Qp_rs = (U.Qp_rs + V.Qp_rs) / 2 + W.p_rs * w.length;
				W.Rs_pq = (U.Rs_pq + V.Rs_pq) / 2 + W.pq_s * w.length;
				W.Sr_pq = (U.Sr_pq + V.Sr_pq) / 2 + W.pq_r * w.length;

				annot.quartetCnt = 0;
				annot.sumLengthA = 0;
				annot.sumLengthB = 0;
				annot.sumLengthC = 0;
				annot.sumLengthD = 0;
				annot.sumInternalLength = 0;
			}

			{
				score_t wp = up + vp, wq = uq + vq, wr = ur + vr, ws = us + vs;

				if (wp < 1e-8) wp = 0;
				if (wq < 1e-8) wq = 0;
				if (wr < 1e-8) wr = 0;
				if (ws < 1e-8) ws = 0;

				if (wp == 0 || wq == 0) W.Pq = W.Qp = W.pqX = W.pq = 0;
				if (wr == 0 || ws == 0) W.Rs = W.Sr = W.rsX = W.rs = 0;

				if (wp == 0 || wq == 0 || wr == 0) W.Sr_pq = W.Qp_r = W.Pq_r = W.pqXr = W.R_pq = W.pq_r = 0;
				if (wp == 0 || wq == 0 || ws == 0) W.Rs_pq = W.Qp_s = W.Pq_s = W.pqXs = W.S_pq = W.pq_s = 0;
				if (wp == 0 || wr == 0 || ws == 0) W.Qp_rs = W.Sr_p = W.Rs_p = W.pXrs = W.P_rs = W.p_rs = 0;
				if (wq == 0 || wr == 0 || ws == 0) W.Pq_rs = W.Sr_q = W.Rs_q = W.qXrs = W.Q_rs = W.q_rs = 0;

				if (wp == 0 || wq == 0 || wr == 0 || ws == 0) annot.sumLengthA = annot.sumLengthB = annot.sumLengthC = annot.sumLengthD = annot.sumInternalLength = annot.quartetCnt = W.pq_rs = 0;
			}
		}

		CustomizedAnnotation extended(Node& w){
			CustomizedAnnotation old = w.annot;

			Node& u = nodes[w.small];
			Node& v = nodes[w.large];

			if (w.dup == false){
				w.aa = u.aa + v.aa;
				w.bb = u.bb + v.bb;
				w.cc = u.cc + v.cc;
				w.dd = u.dd + v.dd;
				
				w.A = u.A + v.A + w.aa * w.length;
				w.B = u.B + v.B + w.bb * w.length;
				w.C = u.C + v.C + w.cc * w.length;
				w.D = u.D + v.D + w.dd * w.length;
			}
			else {
				w.aa = (u.aa + v.aa) / 2;
				w.bb = (u.bb + v.bb) / 2;
				w.cc = (u.cc + v.cc) / 2;
				w.dd = (u.dd + v.dd) / 2;
				
				w.A = (u.A + v.A) / 2 + w.aa * w.length;
				w.B = (u.B + v.B) / 2 + w.bb * w.length;
				w.C = (u.C + v.C) / 2 + w.cc * w.length;
				w.D = (u.D + v.D) / 2 + w.dd * w.length;
			}

			extended<0>(w, u, v);
			extended<1>(w, u, v);
			extended<2>(w, u, v);

			return w.annot - old;
		}
		
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
		CustomizedAnnotation annot;
		
		Partition(const TripartitionInitializer &init, int p): leafParent(init.leafParent[p]), nodes(init.nodes[p].size()), color(init.leafParent[p].size(), -1){
			for (int i = 0; i < nodes.size(); i++){
				nodes[i].up = init.nodes[p][i].up;
				nodes[i].small = init.nodes[p][i].small;
				nodes[i].large = init.nodes[p][i].large;
				nodes[i].dup = init.nodes[p][i].dup;
				nodes[i].length = init.nodes[p][i].length;
			}
		}
		
		void update(int x, int i){
			int y = color[i];
			if (x == y) return;
			if (y != -1) cnt[y]--;
			if (x != -1) cnt[x]++;
			for (int u: leafParent[i]){
				int oldu = u;
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

				u = oldu;
				if (nodes[u].dup == false){
					if (y != -1) ((y == 0) ? nodes[u].aa : (y == 1) ? nodes[u].bb : (y == 2) ? nodes[u].cc : nodes[u].dd)--;
					if (x != -1) ((x == 0) ? nodes[u].aa : (x == 1) ? nodes[u].bb : (x == 2) ? nodes[u].cc : nodes[u].dd)++;
					if (y != -1) ((y == 0) ? nodes[u].A : (y == 1) ? nodes[u].B : (y == 2) ? nodes[u].C : nodes[u].D) -= nodes[u].length;
					if (x != -1) ((x == 0) ? nodes[u].A : (x == 1) ? nodes[u].B : (x == 2) ? nodes[u].C : nodes[u].D) += nodes[u].length;
					w = nodes[u].up;
					while (w != -1){
						CustomizedAnnotation diff = extended(nodes[w]);
						annot += diff;
						w = nodes[w].up;
					}
				}
			}
			color[i] = x;
		}
		
		array<double, 3> score(){
			//LOG << cnt[0] << " " << cnt[1] << " " << cnt[2] << " " << cnt[3] << "; " << totalScore1 << " " << totalScore2 << " " << totalScore3 << "\n";
			double prod = cnt[0] * cnt[1] * cnt[2] * cnt[3];
			return {totalScore1 / prod, totalScore2 / prod, totalScore3 / prod};
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



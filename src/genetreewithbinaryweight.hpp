#define OBJECTIVE_VERSION "2"

/* CHANGE LOG
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
		double sumInternalLength = 0;
		double sumLengthA = 0;
		double sumLengthB = 0;
		double sumLengthC = 0;
		double sumLengthD = 0;

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

	CustomizedAnnotation operator+= (const CustomizedAnnotation& o){
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
		double length = 0;
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
			score_t a = 0, b = 0, c = 0, d = 0;
			score_t ab_c = 0, ab_d = 0, a_cd = 0, b_cd = 0;
			score_t ac_b = 0, a_bd = 0, ac_d = 0, bd_c = 0;
			score_t a_bc = 0, ad_b = 0, ad_c = 0, bc_d = 0;
			score_t abx = 0, acx = 0, adx = 0, bcx = 0, bdx = 0, cdx = 0;
			score_t aby = 0, acy = 0, ady = 0, bcy = 0, bdy = 0, cdy = 0;
			
			score_t ab_cd = 0, ac_bd = 0, ad_bc = 0;
			
			int up = -1, small = -1, large = -1; // -1 for dummy!
			bool isGhostBranch = 1;
			double length = 0;

			double A = 0, B = 0, C = 0, D = 0;
			double Abx = 0, Acx = 0, Adx = 0, Bax = 0, Bcx = 0, Bdx = 0, Cax = 0, Cbx = 0, Cdx = 0, Dax = 0, Dbx = 0, Dcx = 0;
			double Aby = 0, Acy = 0, Ady = 0, Bay = 0, Bcy = 0, Bdy = 0, Cay = 0, Cby = 0, Cdy = 0, Day = 0, Dby = 0, Dcy = 0;
			double Ab_c = 0, Ab_d = 0, Ac_b = 0, Ac_d = 0, Ad_b = 0, Ad_c = 0, A_bc = 0, A_bd = 0, A_cd = 0;
			double Ba_c = 0, Ba_d = 0, Bc_a = 0, Bc_d = 0, Bd_a = 0, Bd_c = 0, B_ac = 0, B_ad = 0, B_cd = 0;
			double Cb_a = 0, Cb_d = 0, Ca_b = 0, Ca_d = 0, Cd_b = 0, Cd_a = 0, C_ab = 0, C_bd = 0, C_ad = 0;
			double Db_c = 0, Db_a = 0, Dc_b = 0, Dc_a = 0, Da_b = 0, Da_c = 0, D_bc = 0, D_ab = 0, D_ac = 0;
			CustomizedAnnotation annot;
		};
		
		array<score_t, 3> normal(Node& w){
			Node& u = nodes[w.small];
			Node& v = nodes[w.large];
			
			w.a = u.a + v.a;
			w.b = u.b + v.b;
			w.c = u.c + v.c;
			w.d = u.d + v.d;
			w.abx = u.abx + v.abx + u.a * v.b + v.a * u.b;
			w.acx = u.acx + v.acx + u.a * v.c + v.a * u.c;
			w.adx = u.adx + v.adx + u.a * v.d + v.a * u.d;
			w.bcx = u.bcx + v.bcx + u.b * v.c + v.b * u.c;
			w.bdx = u.bdx + v.bdx + u.b * v.d + v.b * u.d;
			w.cdx = u.cdx + v.cdx + u.c * v.d + v.c * u.d;
			
			if (w.isGhostBranch) {
				w.aby = u.aby + v.aby + u.a * v.b + v.a * u.b;
				w.acy = u.acy + v.acy + u.a * v.c + v.a * u.c;
				w.ady = u.ady + v.ady + u.a * v.d + v.a * u.d;
				w.bcy = u.bcy + v.bcy + u.b * v.c + v.b * u.c;
				w.bdy = u.bdy + v.bdy + u.b * v.d + v.b * u.d;
				w.cdy = u.cdy + v.cdy + u.c * v.d + v.c * u.d;
			}
			else w.aby = w.acy = w.ady = w.bcy = w.bdy = w.cdy = 0;
			
			w.ab_c = u.ab_c + v.ab_c + u.c * (v.abx - v.aby) + (u.abx - u.aby) * v.c;
			w.ab_d = u.ab_d + v.ab_d + u.d * (v.abx - v.aby) + (u.abx - u.aby) * v.d;
			w.a_cd = u.a_cd + v.a_cd + u.a * (v.cdx - v.cdy) + (u.cdx - u.cdy) * v.a;
			w.b_cd = u.b_cd + v.b_cd + u.b * (v.cdx - v.cdy) + (u.cdx - u.cdy) * v.b;
			
			w.ac_b = u.ac_b + v.ac_b + u.b * (v.acx - v.acy) + (u.acx - u.acy) * v.b;
			w.ac_d = u.ac_d + v.ac_d + u.d * (v.acx - v.acy) + (u.acx - u.acy) * v.d;
			w.a_bd = u.a_bd + v.a_bd + u.a * (v.bdx - v.bdy) + (u.bdx - u.bdy) * v.a;
			w.bd_c = u.bd_c + v.bd_c + u.c * (v.bdx - v.bdy) + (u.bdx - u.bdy) * v.c;
			
			w.ad_b = u.ad_b + v.ad_b + u.b * (v.adx - v.ady) + (u.adx - u.ady) * v.b;
			w.ad_c = u.ad_c + v.ad_c + u.c * (v.adx - v.ady) + (u.adx - u.ady) * v.c;
			w.a_bc = u.a_bc + v.a_bc + u.a * (v.bcx - v.bcy) + (u.bcx - u.bcy) * v.a;
			w.bc_d = u.bc_d + v.bc_d + u.d * (v.bcx - v.bcy) + (u.bcx - u.bcy) * v.d;
			
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
					+ u.bd_c * v.a + v.bd_c * u.a;
			w.ad_bc = u.adx * v.bcx - u.ady * v.bcy + v.adx * u.bcx - v.ady * u.bcy
					+ u.a_bc * v.d + v.a_bc * u.d
					+ u.ad_b * v.c + v.ad_b * u.c
					+ u.ad_c * v.b + v.ad_c * u.b
					+ u.bc_d * v.a + v.bc_d * u.a;

			return {w.ab_cd - ab_cd, w.ac_bd - ac_bd, w.ad_bc - ad_bc};
		}

		CustomizedAnnotation extended(Node& w){
			Node& u = nodes[w.small];
			Node& v = nodes[w.large];
			
			array<score_t, 3> normalArray = normal(w);
			CustomizedAnnotation old = w.annot;
			old.ab_cd.quartetCnt = -normalArray[0];
			old.ac_bd.quartetCnt = -normalArray[1];
			old.ad_bc.quartetCnt = -normalArray[2];

			w.A = u.A + v.A + w.a * w.length;
			w.B = u.B + v.B + w.b * w.length;
			w.C = u.C + v.C + w.c * w.length;
			w.D = u.D + v.D + w.d * w.length;
			
			w.Abx = u.Abx + v.Abx + u.A * v.b + v.A * u.b;
			w.Acx = u.Acx + v.Acx + u.A * v.c + v.A * u.c;
			w.Adx = u.Adx + v.Adx + u.A * v.d + v.A * u.d;
			w.Bax = u.Bax + v.Bax + u.B * v.a + v.B * u.a;
			w.Bcx = u.Bcx + v.Bcx + u.B * v.c + v.B * u.c;
			w.Bdx = u.Bdx + v.Bdx + u.B * v.d + v.B * u.d;
			w.Cax = u.Cax + v.Cax + u.C * v.a + v.C * u.a;
			w.Cbx = u.Cbx + v.Cbx + u.C * v.b + v.C * u.b;
			w.Cdx = u.Cdx + v.Cdx + u.C * v.d + v.C * u.d;
			w.Dax = u.Dax + v.Dax + u.D * v.a + v.D * u.a;
			w.Dbx = u.Dbx + v.Dbx + u.D * v.b + v.D * u.b;
			w.Dcx = u.Dcx + v.Dcx + u.D * v.c + v.D * u.c;
			
			if (w.isGhostBranch) {
				w.Aby = u.Aby + v.Aby + u.A * v.b + v.A * u.b;
				w.Acy = u.Acy + v.Acy + u.A * v.c + v.A * u.c;
				w.Ady = u.Ady + v.Ady + u.A * v.d + v.A * u.d;
				w.Bay = u.Bay + v.Bay + u.B * v.a + v.B * u.a;
				w.Bcy = u.Bcy + v.Bcy + u.B * v.c + v.B * u.c;
				w.Bdy = u.Bdy + v.Bdy + u.B * v.d + v.B * u.d;
				w.Cay = u.Cay + v.Cay + u.C * v.a + v.C * u.a;
				w.Cby = u.Cby + v.Cby + u.C * v.b + v.C * u.b;
				w.Cdy = u.Cdy + v.Cdy + u.C * v.d + v.C * u.d;
				w.Day = u.Day + v.Day + u.D * v.a + v.D * u.a;
				w.Dby = u.Dby + v.Dby + u.D * v.b + v.D * u.b;
				w.Dcy = u.Dcy + v.Dcy + u.D * v.c + v.D * u.c;
			}
			else w.Aby = w.Acy = w.Ady = w.Bay = w.Bcy = w.Bdy = w.Cay = w.Cby = w.Cdy = w.Day = w.Dby = w.Dcy = 0;

			w.Ab_c = u.Ab_c + v.Ab_c + (u.Abx - u.Aby) * v.c + (v.Abx - v.Aby) * u.c;
			w.Ab_d = u.Ab_d + v.Ab_d + (u.Abx - u.Aby) * v.d + (v.Abx - v.Aby) * u.d;
			w.Ac_b = u.Ac_b + v.Ac_b + (u.Acx - u.Acy) * v.b + (v.Acx - v.Acy) * u.b;
			w.Ac_d = u.Ac_d + v.Ac_d + (u.Acx - u.Acy) * v.d + (v.Acx - v.Acy) * u.d;
			w.Ad_b = u.Ad_b + v.Ad_b + (u.Adx - u.Ady) * v.b + (v.Adx - v.Ady) * u.b;
			w.Ad_c = u.Ad_c + v.Ad_c + (u.Adx - u.Ady) * v.c + (v.Adx - v.Ady) * u.c;
			w.Ba_c = u.Ba_c + v.Ba_c + (u.Bax - u.Bay) * v.c + (v.Bax - v.Bay) * u.c;
			w.Ba_d = u.Ba_d + v.Ba_d + (u.Bax - u.Bay) * v.d + (v.Bax - v.Bay) * u.d;
			w.Bc_a = u.Bc_a + v.Bc_a + (u.Bcx - u.Bcy) * v.a + (v.Bcx - v.Bcy) * u.a;
			w.Bc_d = u.Bc_d + v.Bc_d + (u.Bcx - u.Bcy) * v.d + (v.Bcx - v.Bcy) * u.d;
			w.Bd_a = u.Bd_a + v.Bd_a + (u.Bdx - u.Bdy) * v.a + (v.Bdx - v.Bdy) * u.a;
			w.Bd_c = u.Bd_c + v.Bd_c + (u.Bdx - u.Bdy) * v.c + (v.Bdx - v.Bdy) * u.c;
			w.Ca_b = u.Ca_b + v.Ca_b + (u.Cax - u.Cay) * v.b + (v.Cax - v.Cay) * u.b;
			w.Ca_d = u.Ca_d + v.Ca_d + (u.Cax - u.Cay) * v.d + (v.Cax - v.Cay) * u.d;
			w.Cb_a = u.Cb_a + v.Cb_a + (u.Cbx - u.Cby) * v.a + (v.Cbx - v.Cby) * u.a;
			w.Cb_d = u.Cb_d + v.Cb_d + (u.Cbx - u.Cby) * v.d + (v.Cbx - v.Cby) * u.d;
			w.Cd_a = u.Cd_a + v.Cd_a + (u.Cdx - u.Cdy) * v.a + (v.Cdx - v.Cdy) * u.a;
			w.Cd_b = u.Cd_b + v.Cd_b + (u.Cdx - u.Cdy) * v.b + (v.Cdx - v.Cdy) * u.b;
			w.Da_b = u.Da_b + v.Da_b + (u.Dax - u.Day) * v.b + (v.Dax - v.Day) * u.b;
			w.Da_c = u.Da_c + v.Da_c + (u.Dax - u.Day) * v.c + (v.Dax - v.Day) * u.c;
			w.Db_a = u.Db_a + v.Db_a + (u.Dbx - u.Dby) * v.a + (v.Dbx - v.Dby) * u.a;
			w.Db_c = u.Db_c + v.Db_c + (u.Dbx - u.Dby) * v.c + (v.Dbx - v.Dby) * u.c;
			w.Dc_a = u.Dc_a + v.Dc_a + (u.Dcx - u.Dcy) * v.a + (v.Dcx - v.Dcy) * u.a;
			w.Dc_b = u.Dc_b + v.Dc_b + (u.Dcx - u.Dcy) * v.b + (v.Dcx - v.Dcy) * u.b;
			
			w.A_bc = u.A_bc + v.A_bc + u.A * (v.bcx - v.bcy) + v.A * (u.bcx - u.bcy) + w.a_bc * w.length;
			w.A_bd = u.A_bd + v.A_bd + u.A * (v.bdx - v.bdy) + v.A * (u.bdx - u.bdy) + w.a_bd * w.length;
			w.A_cd = u.A_cd + v.A_cd + u.A * (v.cdx - v.cdy) + v.A * (u.cdx - u.cdy) + w.a_cd * w.length;
			w.B_ac = u.B_ac + v.B_ac + u.B * (v.acx - v.acy) + v.B * (u.acx - u.acy) + w.ac_b * w.length;
			w.B_ad = u.B_ad + v.B_ad + u.B * (v.adx - v.ady) + v.B * (u.adx - u.ady) + w.ad_b * w.length;
			w.B_cd = u.B_cd + v.B_cd + u.B * (v.cdx - v.cdy) + v.B * (u.cdx - u.cdy) + w.b_cd * w.length;
			w.C_ab = u.C_ab + v.C_ab + u.C * (v.abx - v.aby) + v.C * (u.abx - u.aby) + w.ab_c * w.length;
			w.C_ad = u.C_ad + v.C_ad + u.C * (v.adx - v.ady) + v.C * (u.adx - u.ady) + w.ad_c * w.length;
			w.C_bd = u.C_bd + v.C_bd + u.C * (v.bdx - v.bdy) + v.C * (u.bdx - u.bdy) + w.bd_c * w.length;
			w.D_ab = u.D_ab + v.D_ab + u.D * (v.abx - v.aby) + v.D * (u.abx - u.aby) + w.ab_d * w.length;
			w.D_ac = u.D_ac + v.D_ac + u.D * (v.acx - v.acy) + v.D * (u.acx - u.acy) + w.ac_d * w.length;
			w.D_bc = u.D_bc + v.D_bc + u.D * (v.bcx - v.bcy) + v.D * (u.bcx - u.bcy) + w.bc_d * w.length;
			
			w.annot.ab_cd.sumLengthA = u.Ab_c * v.d + v.Ab_c * u.d + u.Ab_d * v.c + v.Ab_d * u.c + u.A_cd * v.b + v.A_cd * u.b + u.A * v.b_cd + v.A * u.b_cd + u.Abx * v.cdx - u.Aby * v.cdy + v.Abx * u.cdx - v.Aby * u.cdy;
			w.annot.ab_cd.sumLengthB = u.Ba_c * v.d + v.Ba_c * u.d + u.Ba_d * v.c + v.Ba_d * u.c + u.B_cd * v.a + v.B_cd * u.a + u.B * v.a_cd + v.B * u.a_cd + u.Bax * v.cdx - u.Bay * v.cdy + v.Bax * u.cdx - v.Bay * u.cdy;
			w.annot.ab_cd.sumLengthC = u.Cd_a * v.b + v.Cd_a * u.b + u.Cd_a * v.b + v.Cd_a * u.b + u.C_ab * v.d + v.C_ab * u.d + u.C * v.ab_d + v.C * u.ab_d + u.Cdx * v.abx - u.Cdy * v.aby + v.Cdx * u.abx - v.Cdy * u.aby;
			w.annot.ab_cd.sumLengthD = u.Dc_a * v.b + v.Dc_a * u.b + u.Dc_a * v.b + v.Dc_a * u.b + u.D_ab * v.c + v.D_ab * u.c + u.D * v.ab_c + v.D * u.ab_c + u.Dcx * v.abx - u.Dcy * v.aby + v.Dcx * u.abx - v.Dcy * u.aby;

			w.annot.ac_bd.sumLengthA = u.Ac_b * v.d + v.Ac_b * u.d + u.Ac_d * v.b + v.Ac_d * u.b + u.A_bd * v.c + v.A_bd * u.c + u.A * v.bd_c + v.A * u.bd_c + u.Acx * v.bdx - u.Acy * v.bdy + v.Acx * u.bdx - v.Acy * u.bdy;
            w.annot.ac_bd.sumLengthC = u.Ca_b * v.d + v.Ca_b * u.d + u.Ca_d * v.b + v.Ca_d * u.b + u.C_bd * v.a + v.C_bd * u.a + u.C * v.a_bd + v.C * u.a_bd + u.Cax * v.bdx - u.Cay * v.bdy + v.Cax * u.bdx - v.Cay * u.bdy;
            w.annot.ac_bd.sumLengthB = u.Bd_a * v.c + v.Bd_a * u.c + u.Bd_a * v.c + v.Bd_a * u.c + u.B_ac * v.d + v.B_ac * u.d + u.B * v.ac_d + v.B * u.ac_d + u.Bdx * v.acx - u.Bdy * v.acy + v.Bdx * u.acx - v.Bdy * u.acy;
            w.annot.ac_bd.sumLengthD = u.Db_a * v.c + v.Db_a * u.c + u.Db_a * v.c + v.Db_a * u.c + u.D_ac * v.b + v.D_ac * u.b + u.D * v.ac_b + v.D * u.ac_b + u.Dbx * v.acx - u.Dby * v.acy + v.Dbx * u.acx - v.Dby * u.acy;

			w.annot.ad_bc.sumLengthA = u.Ad_c * v.b + v.Ad_c * u.b + u.Ad_b * v.c + v.Ad_b * u.c + u.A_bc * v.d + v.A_bc * u.d + u.A * v.bc_d + v.A * u.bc_d + u.Adx * v.bcx - u.Ady * v.bcy + v.Adx * u.bcx - v.Ady * u.bcy;
            w.annot.ad_bc.sumLengthD = u.Da_c * v.b + v.Da_c * u.b + u.Da_b * v.c + v.Da_b * u.c + u.D_bc * v.a + v.D_bc * u.a + u.D * v.a_bc + v.D * u.a_bc + u.Dax * v.bcx - u.Day * v.bcy + v.Dax * u.bcx - v.Day * u.bcy;
            w.annot.ad_bc.sumLengthC = u.Cb_a * v.d + v.Cb_a * u.d + u.Cb_a * v.d + v.Cb_a * u.d + u.C_ad * v.b + v.C_ad * u.b + u.C * v.ad_b + v.C * u.ad_b + u.Cbx * v.adx - u.Cby * v.ady + v.Cbx * u.adx - v.Cby * u.ady;
            w.annot.ad_bc.sumLengthB = u.Bc_a * v.d + v.Bc_a * u.d + u.Bc_a * v.d + v.Bc_a * u.d + u.B_ad * v.c + v.B_ad * u.c + u.B * v.ad_c + v.B * u.ad_c + u.Bcx * v.adx - u.Bcy * v.ady + v.Bcx * u.adx - v.Bcy * u.ady;

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
				double dt = t;
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

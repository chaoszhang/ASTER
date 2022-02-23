#include<vector>
#include<array>
#include<thread>

#define SUPPORT

using namespace std;

struct TripartitionInitializer{
	struct Node{
		int up = -1, small = -1, large = -1;
		score_t weight = 1, length = 1;
	};
	
	vector<vector<Node> > nodes;
	vector<vector<vector<int> > > leafParent;
};

struct Tripartition{
	struct Partition{
		struct Node{
			score_t x = 0, x_x = 0;
			score_t y = 0, y_y = 0;
			score_t z = 0, z_z = 0;
			score_t x2a = 0, x2a_x = 0; //const static score_t x2a_x2a = 1;
			score_t y2a = 0, y2a_y = 0; //const static score_t y2a_y2a = 1;
			score_t z2a = 0, z2a_z = 0; //const static score_t z2a_z2a = 1;
			score_t xya = 0, xya_x = 0, xya_y = 0; //const static score_t xya_xya = 1;
			score_t xza = 0, xza_x = 0, xza_z = 0; //const static score_t xza_xza = 1;
			score_t yza = 0, yza_y = 0, yza_z = 0; //const static score_t yza_yza = 1;
			score_t x2b = 0, x2b_x = 0, x2b_x2b = 0;
			score_t y2b = 0, y2b_y = 0, y2b_y2b = 0;
			score_t z2b = 0, z2b_z = 0, z2b_z2b = 0;
			score_t xyb = 0, xyb_x = 0, xyb_y = 0, xyb_xyb = 0;
			score_t xzb = 0, xzb_x = 0, xzb_z = 0, xzb_xzb = 0;
			score_t yzb = 0, yzb_y = 0, yzb_z = 0, yzb_yzb = 0;
			score_t tx = 0, tx_y = 0, tx_z = 0, tx_y2a = 0, tx_z2a = 0, tx_y2b = 0, tx_z2b = 0, tx_tx = 0; 
			score_t ty = 0, ty_x = 0, ty_z = 0, ty_x2a = 0, ty_z2a = 0, ty_x2b = 0, ty_z2b = 0, ty_ty = 0; 
			score_t tz = 0, tz_x = 0, tz_y = 0, tz_x2a = 0, tz_y2a = 0, tz_x2b = 0, tz_y2b = 0, tz_tz = 0; 
			score_t q = 0, q_x = 0, q_y = 0, q_z = 0, q_x2a = 0, q_y2a = 0, q_z2a = 0,
					q_x2b = 0, q_y2b = 0, q_z2b = 0, q_tx = 0, q_ty = 0, q_tz = 0; //const static score_t q_q = 1; 
			
			int up = -1, left = -1, right = -1; // -1 for dummy!
			score_t weight = 1, length = 1;
		};
		
		void normal(Node& w){
			//Node& u = nodes[w.left];
			Node& v = nodes[w.right];
			/*
			w.x = u.x + v.x;
			w.y = u.y + v.y;
			w.z = u.z + v.z;
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
			w.tx = u.tx + v.tx + u.y * (v.z2a - v.z2b) + (u.z2a - u.z2b) * v.y + u.z * (v.y2a - v.y2b) + (u.y2a - u.y2b) * v.z;
			w.ty = u.ty + v.ty + u.x * (v.z2a - v.z2b) + (u.z2a - u.z2b) * v.x + u.z * (v.x2a - v.x2b) + (u.x2a - u.x2b) * v.z;
			w.tz = u.tz + v.tz + u.x * (v.y2a - v.y2b) + (u.y2a - u.y2b) * v.x + u.y * (v.x2a - v.x2b) + (u.x2a - u.x2b) * v.y;
			w.q = u.x * v.tx + u.y * v.ty + u.z * v.tz + v.x * u.tx + v.y * u.ty + v.z * u.tz + u.q + v.q
				+ u.x2a * v.yza - u.x2b * v.yzb
				+ u.y2a * v.xza - u.y2b * v.xzb
				+ u.z2a * v.xya - u.z2b * v.xyb;
			*/
			w.x = v.x * w.length; w.x_x = w.length;
			w.y = v.y * w.length; w.y_y = w.length;
			w.z = v.z * w.length; w.z_z = w.length;
			w.x2a = v.x2a; w.x2a_x = v.x; // 1
			w.y2a = v.y2a; w.y2a_y = v.y;
			w.z2a = v.z2a; w.z2a_z = v.z;
			w.xya = v.xya; w.xya_x = v.y; w.xya_y = v.x; // 1
			w.xza = v.xza; w.xza_x = v.z; w.xza_z = v.x;
			w.yza = v.yza; w.yza_y = v.z; w.yza_z = v.y;
			w.x2b = v.x2b * (1 - w.weight); w.x2b_x = v.x * (1 - w.weight); w.x2b_x2b = 1 - w.weight;
			w.y2b = v.y2b * (1 - w.weight); w.y2b_y = v.y * (1 - w.weight); w.y2b_y2b = 1 - w.weight;
			w.z2b = v.z2b * (1 - w.weight); w.z2b_z = v.z * (1 - w.weight); w.z2b_z2b = 1 - w.weight;
			w.xyb = v.xyb * (1 - w.weight); w.xyb_x = v.y * (1 - w.weight); w.xyb_y = v.x * (1 - w.weight); w.xyb_xyb = 1 - w.weight;
			w.xzb = v.xzb * (1 - w.weight); w.xzb_x = v.z * (1 - w.weight); w.xzb_z = v.x * (1 - w.weight); w.xzb_xzb = 1 - w.weight;
			w.yzb = v.yzb * (1 - w.weight); w.yzb_y = v.z * (1 - w.weight); w.yzb_z = v.y * (1 - w.weight); w.yzb_yzb = 1 - w.weight;
			w.tx = v.tx * w.length; w.tx_y = (v.z2a - v.z2b) * w.length; w.tx_z = (v.y2a - v.y2b) * w.length; w.tx_y2a = v.z * w.length; w.tx_z2a = v.y * w.length; w.tx_y2b = -v.z * w.length; w.tx_z2b = -v.y * w.length; w.tx_tx = w.length;
			w.ty = v.ty * w.length; w.ty_x = (v.z2a - v.z2b) * w.length; w.ty_z = (v.y2a - v.y2b) * w.length; w.ty_x2a = v.z * w.length; w.ty_z2a = v.x * w.length; w.ty_x2a = -v.z * w.length; w.ty_z2a = -v.x * w.length; w.ty_ty = w.length;
			w.tz = v.tz * w.length; w.tz_x = (v.y2a - v.y2b) * w.length; w.tz_y = (v.x2a - v.x2b) * w.length; w.tz_x2a = v.y * w.length; w.tz_y2a = v.x * w.length; w.tz_x2a = -v.y * w.length; w.tz_y2a = -v.x * w.length; w.tz_tz = w.length;
			
			w.q = v.q; w.q_x = v.tx; w.q_y = v.ty; w.q_z = v.tz; w.q_x2a = v.yza; w.q_y2a = v.xza; w.q_z2a = v.xya;
			w.q_x2b = -v.yzb; w.q_y2b = -v.xzb; w.q_z2b = -v.xyb; w.q_tx = v.x; w.q_ty = v.y; w.q_tz = v.z; // 1
		}
		
		void special(Node& w){
			Node& u = nodes[w.left];
			Node& v = nodes[w.right];
			/*
			w.x = u.x + v.x;
			w.y = u.y + v.y;
			w.z = u.z + v.z;
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
			w.tx = u.tx + v.tx + u.y * (v.z2a - v.z2b) + (u.z2a - u.z2b) * v.y + u.z * (v.y2a - v.y2b) + (u.y2a - u.y2b) * v.z;
			w.ty = u.ty + v.ty + u.x * (v.z2a - v.z2b) + (u.z2a - u.z2b) * v.x + u.z * (v.x2a - v.x2b) + (u.x2a - u.x2b) * v.z;
			w.tz = u.tz + v.tz + u.x * (v.y2a - v.y2b) + (u.y2a - u.y2b) * v.x + u.y * (v.x2a - v.x2b) + (u.x2a - u.x2b) * v.y;
			w.q = u.x * v.tx + u.y * v.ty + u.z * v.tz + v.x * u.tx + v.y * u.ty + v.z * u.tz + u.q + v.q
				+ u.x2a * v.yza - u.x2b * v.yzb
				+ u.y2a * v.xza - u.y2b * v.xzb
				+ u.z2a * v.xya - u.z2b * v.xyb;
			*/
			w.x = v.x*1 + v.x_x*u.x; w.x_x = v.x_x*u.x_x;
			w.y = v.y*1 + v.y_y*u.y; w.y_y = v.y_y*u.y_y;
			w.z = v.z*1 + v.z_z*u.z; w.z_z = v.z_z*u.z_z;
			
			w.x2a = v.x2a*1 + v.x2a_x*u.x + 1*u.x2a; w.x2a_x = v.x2a_x*u.x_x + 1*u.x2a_x; // 1
			w.y2a = v.y2a*1 + v.y2a_y*u.y + 1*u.y2a; w.y2a_y = v.y2a_y*u.y_y + 1*u.y2a_y;
			w.z2a = v.z2a*1 + v.z2a_z*u.z + 1*u.z2a; w.z2a_z = v.z2a_z*u.z_z + 1*u.z2a_z;
			w.xya = v.xya*1 + v.xya_x*u.x + v.xya_y*u.y + 1*u.xya; w.xya_x = v.xya_x*u.x_x + u.xya_x; w.xya_y = v.xya_y*u.y_y + 1*u.xya_y; // 1
			w.xza = v.xza*1 + v.xza_x*u.x + v.xza_z*u.z + 1*u.xza; w.xza_x = v.xza_x*u.x_x + u.xza_x; w.xza_z = v.xza_z*u.z_z + 1*u.xza_z;
			w.yza = v.yza*1 + v.yza_y*u.y + v.yza_z*u.z + 1*u.yza; w.yza_y = v.yza_y*u.y_y + u.yza_y; w.yza_z = v.yza_z*u.z_z + 1*u.yza_z;
			w.x2b = v.x2b*1 + v.x2b_x*u.x + v.x2b_x2b*u.x2b; w.x2b_x = v.x2b_x*u.x_x + v.x2b_x2b*u.x2b_x; w.x2b_x2b = v.x2b_x2b*u.x2b_x2b;
			w.y2b = v.y2b*1 + v.y2b_y*u.y + v.y2b_y2b*u.y2b; w.y2b_y = v.y2b_y*u.y_y + v.y2b_y2b*u.y2b_y; w.y2b_y2b = v.y2b_y2b*u.y2b_y2b;
			w.z2b = v.z2b*1 + v.z2b_z*u.z + v.z2b_z2b*u.z2b; w.z2b_z = v.z2b_z*u.z_z + v.z2b_z2b*u.z2b_z; w.z2b_z2b = v.z2b_z2b*u.z2b_z2b;
			w.xyb = v.xyb*1 + v.xyb_x*u.x + v.xyb_y*u.y + v.xyb_xyb*u.xyb; w.xyb_x = v.xyb_x*u.x_x + v.xyb_xyb*u.xyb_x; w.xyb_y = v.xyb_y*u.y_y + v.xyb_xyb*u.xyb_y; w.xyb_xyb = v.xyb_xyb*u.xyb_xyb;
			w.xzb = v.xzb*1 + v.xzb_x*u.x + v.xzb_z*u.z + v.xzb_xzb*u.xzb; w.xzb_x = v.xzb_x*u.x_x + v.xzb_xzb*u.xzb_x; w.xzb_z = v.xzb_z*u.z_z + v.xzb_xzb*u.xzb_z; w.xzb_xzb = v.xzb_xzb*u.xzb_xzb;
			w.yzb = v.yzb*1 + v.yzb_y*u.y + v.yzb_z*u.z + v.yzb_yzb*u.yzb; w.yzb_y = v.yzb_y*u.y_y + v.yzb_yzb*u.yzb_y; w.yzb_z = v.yzb_z*u.z_z + v.yzb_yzb*u.yzb_z; w.yzb_yzb = v.yzb_yzb*u.yzb_yzb;
			
			w.tx = v.tx*1 + v.tx_y*u.y + v.tx_z*u.z + v.tx_y2a*u.y2a + v.tx_z2a*u.z2a + v.tx_tx*u.tx;
			w.tx_y = v.tx_y*u.y_y + v.tx_y2a*u.y2a_y + v.tx_y2b*u.y2b_y + v.tx_tx*u.tx_y; w.tx_z = v.tx_z*u.z_z + v.tx_z2a*u.z2a_z + v.tx_z2b*u.z2b_z + v.tx_tx*u.tx_z;
			w.tx_y2a = v.tx_y2a*1 + v.tx_tx*u.tx_y2a; w.tx_z2a = v.tx_z2a*1 + v.tx_tx*u.tx_z2a; w.tx_y2b = v.tx_y2b*u.y2b_y2b + v.tx_tx*u.tx_y2b; w.tx_z2b = v.tx_z2b*u.z2b_z2b + v.tx_tx*u.tx_z2b;
			w.tx_tx = v.tx_tx*u.tx_tx;
			
			w.ty = v.ty*1 + v.ty_x*u.x + v.ty_z*u.z + v.ty_x2a*u.x2a + v.ty_z2a*u.z2a + v.ty_ty*u.ty;
			w.ty_x = v.ty_x*u.x_x + v.ty_x2a*u.x2a_x + v.ty_x2b*u.x2b_x + v.ty_ty*u.ty_x; w.ty_z = v.ty_z*u.z_z + v.ty_z2a*u.z2a_z + v.ty_z2b*u.z2b_z + v.ty_ty*u.ty_z;
			w.ty_x2a = v.ty_x2a*1 + v.ty_ty*u.ty_x2a; w.ty_z2a = v.ty_z2a*1 + v.ty_ty*u.ty_z2a; w.ty_x2b = v.ty_x2b*u.x2b_x2b + v.ty_ty*u.ty_x2b; w.ty_z2b = v.ty_z2b*u.z2b_z2b + v.ty_ty*u.ty_z2b;
			w.ty_ty = v.ty_ty*u.ty_ty;
			
			w.tz = v.tz*1 + v.tz_x*u.x + v.tz_y*u.y + v.tz_x2a*u.x2a + v.tz_y2a*u.y2a + v.tz_tz*u.tz;
			w.tz_x = v.tz_x*u.x_x + v.tz_x2a*u.x2a_x + v.tz_x2b*u.x2b_x + v.tz_tz*u.tz_x; w.tz_y = v.tz_y*u.y_y + v.tz_y2a*u.y2a_y + v.tz_y2b*u.y2b_y + v.tz_tz*u.tz_y;
			w.tz_x2a = v.tz_x2a*1 + v.tz_tz*u.tz_x2a; w.tz_y2a = v.tz_y2a*1 + v.tz_tz*u.tz_y2a; w.tz_x2b = v.tz_x2b*u.x2b_x2b + v.tz_tz*u.tz_x2b; w.tz_y2b = v.tz_y2b*u.y2b_y2b + v.tz_tz*u.tz_y2b;
			w.tz_tz = v.tz_tz*u.tz_tz;
			
			w.q = v.q*1 + v.q_x*u.x + v.q_y*u.y + v.q_z*u.z + v.q_x2a*u.x2a + v.q_y2a*u.y2a + v.q_z2a*u.z2a + v.q_x2b*u.x2b + v.q_y2b*u.y2b + v.q_z2b*u.z2b + v.q_tx*u.tx + v.q_ty*u.ty + v.q_tz*u.tz + 1*u.q;
			w.q_x = v.q_x*u.x_x + v.q_x2a*u.x2a_x + v.q_x2b*u.x2b_x + v.q_ty*u.ty_x + v.q_tz*u.tz_x + 1*u.q_x;
			w.q_y = v.q_y*u.y_y + v.q_y2a*u.y2a_y + v.q_y2b*u.y2b_y + v.q_tx*u.tx_y + v.q_tz*u.tz_y + 1*u.q_y; 
			w.q_z = v.q_z*u.z_z + v.q_z2a*u.z2a_z + v.q_z2b*u.z2b_z + v.q_tx*u.tx_z + v.q_ty*u.ty_z + 1*u.q_z; 
			w.q_x2a = v.q_x2a*1 + v.q_ty*u.ty_x2a + v.q_tz*u.tz_x2a + 1*u.q_x2a; w.q_x2b = v.q_x2b*u.x2b_x2b + v.q_ty*u.ty_x2b + v.q_tz*u.tz_x2b + 1*u.q_x2b;
			w.q_y2a = v.q_y2a*1 + v.q_tx*u.tx_y2a + v.q_tz*u.tz_y2a + 1*u.q_y2a; w.q_y2b = v.q_y2b*u.y2b_y2b + v.q_tx*u.tx_y2b + v.q_tz*u.tz_y2b + 1*u.q_y2b;
			w.q_z2a = v.q_z2a*1 + v.q_tx*u.tx_z2a + v.q_ty*u.ty_z2a + 1*u.q_z2a; w.q_z2b = v.q_z2b*u.z2b_z2b + v.q_tx*u.tx_z2b + v.q_ty*u.ty_z2b + 1*u.q_z2b;
			w.q_tx = v.q_tx * u.tx_tx + 1*u.q_tx; w.q_ty = v.q_ty * u.ty_ty + 1*u.q_ty; w.q_tz = v.q_tz * u.tz_tz + 1*u.q_tz;
		}
		
		vector<vector<int> > leafParent;
		score_t totalScore = 0;
		vector<Node> nodes;
		vector<int> color;
		
		void preprocess(int v, const vector<TripartitionInitializer::Node> &inodes, vector<int> &size){
			if (size[v] > 0) return;
			if (inodes[v].small == -1 && inodes[v].large == -1) size[v] = 1;
			else {
				preprocess(inodes[v].small, inodes, size);
				preprocess(inodes[v].large, inodes, size);
				size[v] = size[inodes[v].small] + size[inodes[v].large];
			}
		}
		
		int initialize(int v, int u, const vector<TripartitionInitializer::Node> &inodes, vector<int> &size, vector<int> &leafId){
			const int cur = nodes.size();
			cerr << cur << " " << v << " " << u << " " << inodes[v].large << " " << inodes[v].small << endl;
			nodes.emplace_back();
			if (leafId[v] != -1) {
				cerr << "case 1\n";
				leafParent[leafId[v]].push_back(cur);
				nodes[cur].length = inodes[v].length;
				return cur;
			}
			if (((size[inodes[v].large] > size[inodes[v].small]) ? inodes[v].large : inodes[v].small) == u){
				cerr << "case 2\n";
				nodes[cur].left = -1;
				nodes[cur].right = initialize((size[inodes[v].large] > size[inodes[v].small]) ? inodes[v].small : inodes[v].large, -1, inodes, size, leafId);
				nodes[nodes[cur].right].up = cur;
				nodes[cur].weight = inodes[v].weight;
				nodes[cur].length = inodes[v].length;
				cerr << "case 2' " << cur << "\n";
				return cur;
			}
			cerr << "case 3\n";
			int w = (size[inodes[v].large] > size[inodes[v].small]) ? inodes[v].large : inodes[v].small, su = (u == -1) ? 0 : size[u];
			while (((size[inodes[w].large] > size[inodes[w].small]) ? inodes[w].large : inodes[w].small) != u && 3 * (size[w] - su) > 2 * (size[v] - su)){
				w = (size[inodes[w].large] > size[inodes[w].small]) ? inodes[w].large : inodes[w].small;
			}
			nodes[cur].left = initialize(w, u, inodes, size, leafId);
			cerr << "case 3'\n";
			nodes[cur].right = initialize(v, w, inodes, size, leafId);
			cerr << "case 3'' " << nodes[cur].left << " " << nodes[cur].right << endl;
			nodes[nodes[cur].left].up = cur;
			nodes[nodes[cur].right].up = cur;
			cerr << "case 3'''\n";
			return cur;
		}
		
		Partition(const TripartitionInitializer &init, int p): leafParent(init.leafParent[p].size()), color(init.leafParent[p].size(), -1){
			vector<int> size(init.nodes[p].size()), leafId(init.nodes[p].size(), -1);
			
			for (int i = 0; i < init.nodes[p].size(); i++) preprocess(i, init.nodes[p], size);
			for (int i = 0; i < init.leafParent[p].size(); i++){
				for (int v: init.leafParent[p][i]) leafId[v] = i;
			}
			for (int i = 0; i < init.nodes[p].size(); i++){
				if (init.nodes[p][i].up == -1) initialize(i, -1, init.nodes[p], size, leafId);
			}
		}
		
		void update(int x, int i){
			int y = color[i], q;
			if (x == y) return;
			for (int u: leafParent[i]){
				if (y != -1) ((y == 0) ? nodes[u].z : (y == 1) ? nodes[u].x : nodes[u].y) -= nodes[u].length;
				if (x != -1) ((x == 0) ? nodes[u].z : (x == 1) ? nodes[u].x : nodes[u].y) += nodes[u].length;
				int w = nodes[u].up;
				while (w != -1){
					q = -nodes[w].q;
					if (nodes[w].left == -1) normal(nodes[w]);
					else special(nodes[w]);
					q += nodes[w].q;
					u = w;
					w = nodes[u].up;
				}
				totalScore += q;
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
		vector<score_t> score1, score2, score3;
		score_t totalScore1 = 0, totalScore2 = 0, totalScore3 = 0, cnt[4] = {};
		vector<Node> nodes;
		vector<int> color;
		
		Partition(TripartitionInitializer init, int p): leafParent(init.leafParent[p]), nodes(init.nodes[p].size()), 
				score1(init.nodes[p].size()), score2(init.nodes[p].size()), score3(init.nodes[p].size()), color(init.leafParent[p].size(), -1){
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
			if (y != -1) cnt[y]--;
			if (x != -1) cnt[x]++;
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
			}
			color[i] = x;
		}
		
		array<double, 3> score(){
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

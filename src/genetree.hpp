#define OBJECTIVE_VERSION "1"

/* CHANGE LOG
 * 1. using thread pool for better parallelism
 */

#include <vector>
#include <array>
#include <thread>
#include "threadpool.hpp"

using namespace std;

struct TripartitionInitializer{
	struct Node{
		int up = -1, down = -1; // -1 for dummy!
		int maxChild = -1, max2Child = -1, downMost = -1;
		int maxHeight = 0, max2Height = 0, path = 0;
	};
	
	int nTaxa;
	vector<vector<int> > roots;
	vector<vector<vector<pair<int, int> > > > parentChild;
};

struct Tripartition{
	struct Partition{
		struct Node{
			score_t X = 0, Y = 0, Z = 0,
					sX = 0, sY = 0, sZ = 0, sXY = 0, sYZ = 0, sXZ = 0, sX2 = 0, sY2 = 0, sZ2 = 0,
					sX2Y = 0, sX2Z = 0, sY2X = 0, sY2Z = 0, sZ2X = 0, sZ2Y = 0, sX2YZ = 0, sY2XZ = 0, sZ2XY = 0,
					sTx = 0, sTxX = 0, sTy = 0, sTyY = 0, sTz = 0, sTzZ = 0, sQ = 0,
					Tx = 0, Txuy = 0, Txuz = 0,
					Ty = 0, Tyux = 0, Tyuz = 0,
					Tz = 0, Tzux = 0, Tzuy = 0,
					Q = 0, Qux = 0, Quy = 0, Quz = 0, Quxy = 0, Quxz = 0, Quyz = 0, Qux2 = 0, Quy2 = 0, Quz2 = 0; //QuQ === 1
			int up = -1, down = -1; // -1 for dummy!
			
			void testAdd(Node &v){
				const score_t vX2 = v.X*(v.X-1)/2, vY2 = v.Y*(v.Y-1)/2, vZ2 = v.Z*(v.Z-1)/2;
			
				Node &u = *this;
				u.sX += v.X;
				u.sY += v.Y;
				u.sZ += v.Z;
				u.sXY += v.X * v.Y;
				u.sXZ += v.X * v.Z;
				u.sYZ += v.Y * v.Z;
				u.sX2 += vX2;
				u.sY2 += vY2;
				u.sZ2 += vZ2;
				u.sX2Y += vX2 * v.Y;
				u.sX2Z += vX2 * v.Z;
				u.sY2X += vY2 * v.X;
				u.sY2Z += vY2 * v.Z;
				u.sZ2X += vZ2 * v.X;
				u.sZ2Y += vZ2 * v.Y;
				u.sX2YZ += vX2 * v.Y * v.Z;
				u.sY2XZ += vY2 * v.X * v.Z;
				u.sZ2XY += vZ2 * v.X * v.Y;
				u.sTx += v.Tx; 
				u.sTxX += v.Tx * v.X;
				u.sTy += v.Ty;
				u.sTyY += v.Ty * v.Y;
				u.sTz += v.Tz;
				u.sTzZ += v.Tz * v.Z;
				u.sQ += v.Q;
			}
			
			void testRmv(Node &v){
				const score_t vX2 = v.X*(v.X-1)/2, vY2 = v.Y*(v.Y-1)/2, vZ2 = v.Z*(v.Z-1)/2;
			
				Node &u = *this;
				u.sX -= v.X;
				u.sY -= v.Y;
				u.sZ -= v.Z;
				u.sXY -= v.X * v.Y;
				u.sXZ -= v.X * v.Z;
				u.sYZ -= v.Y * v.Z;
				u.sX2 -= vX2;
				u.sY2 -= vY2;
				u.sZ2 -= vZ2;
				u.sX2Y -= vX2 * v.Y;
				u.sX2Z -= vX2 * v.Z;
				u.sY2X -= vY2 * v.X;
				u.sY2Z -= vY2 * v.Z;
				u.sZ2X -= vZ2 * v.X;
				u.sZ2Y -= vZ2 * v.Y;
				u.sX2YZ -= vX2 * v.Y * v.Z;
				u.sY2XZ -= vY2 * v.X * v.Z;
				u.sZ2XY -= vZ2 * v.X * v.Y;
				u.sTx -= v.Tx; 
				u.sTxX -= v.Tx * v.X;
				u.sTy -= v.Ty;
				u.sTyY -= v.Ty * v.Y;
				u.sTz -= v.Tz;
				u.sTzZ -= v.Tz * v.Z;
				u.sQ -= v.Q;
			}
			
			void testAddSpecial(Node &w){
				Node &u = *this;
				score_t usX_c2 = u.sX*(u.sX-1)/2, usY_c2 = u.sY*(u.sY-1)/2, usZ_c2 = u.sZ*(u.sZ-1)/2;
				
				u.X += w.X;
				u.Y += w.Y;
				u.Z += w.Z;
				u.Q += w.Q + w.Qux*u.sX + w.Quy*u.sY + w.Quz*u.sZ + w.Quxy*u.sX*u.sY + w.Quxz*u.sX*u.sZ + w.Quyz*u.sY*u.sZ + w.Qux2*usX_c2 + w.Quy2*usY_c2 + w.Quz2*usZ_c2
					+ 2*w.X*(u.sTx + u.sY*u.sZ2 + u.sY2*u.sZ - u.sZ2Y - u.sY2Z) + 2*w.Y*(u.sTy + u.sX*u.sZ2 + u.sX2*u.sZ - u.sZ2X - u.sX2Z) + 2*w.Z*(u.sTz + u.sX*u.sY2 + u.sX2*u.sY - u.sX2Y - u.sY2X);
				u.Qux += 2*w.Y*u.sZ2 + 2*w.Z*u.sY2 + w.Qux + w.Qux2*u.sX + w.Quxy*u.sY + w.Quxz*u.sZ;
				u.Quy += 2*w.X*u.sZ2 + 2*w.Z*u.sX2 + w.Quy + w.Quxy*u.sX + w.Quy2*u.sY + w.Quyz*u.sZ;
				u.Quz += 2*w.X*u.sY2 + 2*w.Y*u.sX2 + w.Quz + w.Quxz*u.sX + w.Quyz*u.sY + w.Quz2*u.sZ;
				u.Quxy += w.Quxy;
				u.Quxz += w.Quxz;
				u.Quyz += w.Quyz;
				u.Qux2 += 2*w.Y*u.sZ + 2*w.Z*u.sY + w.Qux2;
				u.Quy2 += 2*w.X*u.sZ + 2*w.Z*u.sX + w.Quy2;
				u.Quz2 += 2*w.X*u.sY + 2*w.Y*u.sX + w.Quz2;
				u.Tx += w.Tx + w.Txuy*u.sY + w.Z*usY_c2 + w.Txuz*u.sZ + w.Y*usZ_c2;
				u.Txuy += w.Txuy + w.Z*u.sY;
				u.Txuz += w.Txuz + w.Y*u.sZ;
				u.Ty += w.Ty + w.Tyux*u.sX + w.Z*usX_c2 + w.Tyuz*u.sZ + w.X*usZ_c2;
				u.Tyux += w.Tyux + w.Z*u.sX;
				u.Tyuz += w.Tyuz + w.X*u.sZ;
				u.Tz += w.Tz + w.Tzux*u.sX + w.Y*usX_c2 + w.Tzuy*u.sY + w.X*usY_c2;
				u.Tzux += w.Tzux + w.Y*u.sX;
				u.Tzuy += w.Tzuy + w.X*u.sY;
				
			}
			
			void testRmvSpecial(Node &w){
				Node &u = *this;
				score_t usX_c2 = u.sX*(u.sX-1)/2, usY_c2 = u.sY*(u.sY-1)/2, usZ_c2 = u.sZ*(u.sZ-1)/2;
				
				u.X -= w.X;
				u.Y -= w.Y;
				u.Z -= w.Z;
				u.Q -= w.Q + w.Qux*u.sX + w.Quy*u.sY + w.Quz*u.sZ + w.Quxy*u.sX*u.sY + w.Quxz*u.sX*u.sZ + w.Quyz*u.sY*u.sZ + w.Qux2*usX_c2 + w.Quy2*usY_c2 + w.Quz2*usZ_c2
					+ 2*w.X*(u.sTx + u.sY*u.sZ2 + u.sY2*u.sZ - u.sZ2Y - u.sY2Z) + 2*w.Y*(u.sTy + u.sX*u.sZ2 + u.sX2*u.sZ - u.sZ2X - u.sX2Z) + 2*w.Z*(u.sTz + u.sX*u.sY2 + u.sX2*u.sY - u.sX2Y - u.sY2X);
				u.Qux -= 2*w.Y*u.sZ2 + 2*w.Z*u.sY2 + w.Qux + w.Qux2*u.sX + w.Quxy*u.sY + w.Quxz*u.sZ;
				u.Quy -= 2*w.X*u.sZ2 + 2*w.Z*u.sX2 + w.Quy + w.Quxy*u.sX + w.Quy2*u.sY + w.Quyz*u.sZ;
				u.Quz -= 2*w.X*u.sY2 + 2*w.Y*u.sX2 + w.Quz + w.Quxz*u.sX + w.Quyz*u.sY + w.Quz2*u.sZ;
				u.Quxy -= w.Quxy;
				u.Quxz -= w.Quxz;
				u.Quyz -= w.Quyz;
				u.Qux2 -= 2*w.Y*u.sZ + 2*w.Z*u.sY + w.Qux2;
				u.Quy2 -= 2*w.X*u.sZ + 2*w.Z*u.sX + w.Quy2;
				u.Quz2 -= 2*w.X*u.sY + 2*w.Y*u.sX + w.Quz2;
				u.Tx -= w.Tx + w.Txuy*u.sY + w.Z*usY_c2 + w.Txuz*u.sZ + w.Y*usZ_c2;
				u.Txuy -= w.Txuy + w.Z*u.sY;
				u.Txuz -= w.Txuz + w.Y*u.sZ;
				u.Ty -= w.Ty + w.Tyux*u.sX + w.Z*usX_c2 + w.Tyuz*u.sZ + w.X*usZ_c2;
				u.Tyux -= w.Tyux + w.Z*u.sX;
				u.Tyuz -= w.Tyuz + w.X*u.sZ;
				u.Tz -= w.Tz + w.Tzux*u.sX + w.Y*usX_c2 + w.Tzuy*u.sY + w.X*usY_c2;
				u.Tzux -= w.Tzux + w.Y*u.sX;
				u.Tzuy -= w.Tzuy + w.X*u.sY;
			}
			
			void testCalc(Node &w){
				Node &u = *this;
			/*
			3*sx2yz + 3*sy2xz + 3*sz2xy - 2*sx*sy2z - 2*sx2z*sy - sx2*syz - sxz*sy2 - 2*sx*sz2y - 2*sx2y*sz - sxy*sz2 - 2*sy*sz2x - 2*sy2x*sz + 2*sx*sy*sz2 + 2*sx*sy2*sz + 2*sx2*sy*sz
			+ 2*sy*sz2*ux + 2*sy2*sz*ux - 2*sy2z*ux - 2*sz2y*ux 
			+ 2*sy*sz*ux2 - syz*ux2
			+ 2*sx*sz2*uy + 2*sx2*sz*uy - 2*sx2z*uy - 2*sz2x*uy
			+ 2*sx*sz*uy2 - sxz*uy2
			+ 2*sx*sy2*uz + 2*sx2*sy*uz - 2*sx2y*uz - 2*sy2x*uz
			+ 2*sx*sy*uz2 - sxy*uz2 
			+ sx2*uy*uz + sy2*ux*uz + sz2*ux*uy
			tx = sy*sz2 + sy2*sz - sz2y - sy2z + sy2*uz + sz2*uy + sy*uz2 + sz*uy2 + stx + utx
			ty = sx*sz2 + sx2*sz - sz2x - sx2z + sx2*uz + sz2*ux + sx*uz2 + sz*ux2 + sty + uty
			tz = sx*sy2 + sx2*sy - sy2x - sx2y + sx2*uy + sy2*ux + sx*uy2 + sy*ux2 + stz + utz
			Tx = tx + w.Tx
			Ty = ty + w.Ty
			Tz = tz + w.Tz
			
			ux -> sx + ux
			ux2 -> sx_c2 + sx * ux + ux2
			uxy -> sx * sy + sx * uy + sy * ux + uxy
			utx -> sy*sz2 + sy2*sz - sz2y - sy2z + sy2*uz + sz2*uy + sy*uz2 + sz*uy2 + stx + utx
			 */
				score_t usX_c2 = u.sX*(u.sX-1)/2, usY_c2 = u.sY*(u.sY-1)/2, usZ_c2 = u.sZ*(u.sZ-1)/2;
				
				u.X = u.sX + w.X;
				u.Y = u.sY + w.Y;
				u.Z = u.sZ + w.Z;
				u.Q = u.sQ + w.Q + w.Qux*u.sX + w.Quy*u.sY + w.Quz*u.sZ + w.Quxy*u.sX*u.sY + w.Quxz*u.sX*u.sZ + w.Quyz*u.sY*u.sZ + w.Qux2*usX_c2 + w.Quy2*usY_c2 + w.Quz2*usZ_c2
					+ 2*w.X*(u.sTx + u.sY*u.sZ2 + u.sY2*u.sZ - u.sZ2Y - u.sY2Z) + 2*w.Y*(u.sTy + u.sX*u.sZ2 + u.sX2*u.sZ - u.sZ2X - u.sX2Z) + 2*w.Z*(u.sTz + u.sX*u.sY2 + u.sX2*u.sY - u.sX2Y - u.sY2X)
					+ 3*u.sX2YZ + 3*u.sY2XZ + 3*u.sZ2XY - 2*u.sX*u.sY2Z - 2*u.sX2Z*u.sY - u.sX2*u.sYZ - u.sXZ*u.sY2 - 2*u.sX*u.sZ2Y - 2*u.sX2Y*u.sZ - u.sXY*u.sZ2 - 2*u.sY*u.sZ2X
					- 2*u.sY2X*u.sZ + 2*u.sX*u.sY*u.sZ2 + 2*u.sX*u.sY2*u.sZ + 2*u.sX2*u.sY*u.sZ
					+ 2*u.sTx*u.sX - 2*u.sTxX + 2*u.sTy*u.sY - 2*u.sTyY + 2*u.sTz*u.sZ - 2*u.sTzZ;
				u.Qux = 2*w.Y*u.sZ2 + 2*w.Z*u.sY2 + w.Qux + w.Qux2*u.sX + w.Quxy*u.sY + w.Quxz*u.sZ + 2*u.sY*u.sZ2 + 2*u.sY2*u.sZ - 2*u.sY2Z - 2*u.sZ2Y + 2*u.sTx;
				u.Quy = 2*w.X*u.sZ2 + 2*w.Z*u.sX2 + w.Quy + w.Quxy*u.sX + w.Quy2*u.sY + w.Quyz*u.sZ + 2*u.sX*u.sZ2 + 2*u.sX2*u.sZ - 2*u.sX2Z - 2*u.sZ2X + 2*u.sTy;
				u.Quz = 2*w.X*u.sY2 + 2*w.Y*u.sX2 + w.Quz + w.Quxz*u.sX + w.Quyz*u.sY + w.Quz2*u.sZ + 2*u.sX*u.sY2 + 2*u.sX2*u.sY - 2*u.sX2Y - 2*u.sY2X + 2*u.sTz;
				u.Quxy = w.Quxy + u.sZ2;
				u.Quxz = w.Quxz + u.sY2;
				u.Quyz = w.Quyz + u.sX2;
				u.Qux2 = 2*w.Y*u.sZ + 2*w.Z*u.sY + w.Qux2 + 2*u.sY*u.sZ - u.sYZ;
				u.Quy2 = 2*w.X*u.sZ + 2*w.Z*u.sX + w.Quy2 + 2*u.sX*u.sZ - u.sXZ;
				u.Quz2 = 2*w.X*u.sY + 2*w.Y*u.sX + w.Quz2 + 2*u.sX*u.sY - u.sXY;
				u.Tx = w.Tx + w.Txuy*u.sY + w.Z*usY_c2 + w.Txuz*u.sZ + w.Y*usZ_c2 + u.sTx + u.sY*u.sZ2 + u.sY2*u.sZ - u.sZ2Y - u.sY2Z;
				u.Txuy = w.Txuy + w.Z*u.sY + u.sZ2;
				u.Txuz = w.Txuz + w.Y*u.sZ + u.sY2;
				u.Ty = w.Ty + w.Tyux*u.sX + w.Z*usX_c2 + w.Tyuz*u.sZ + w.X*usZ_c2 + u.sTy + u.sX*u.sZ2 + u.sX2*u.sZ - u.sZ2X - u.sX2Z;
				u.Tyux = w.Tyux + w.Z*u.sX + u.sZ2;
				u.Tyuz = w.Tyuz + w.X*u.sZ + u.sX2;
				u.Tz = w.Tz + w.Tzux*u.sX + w.Y*usX_c2 + w.Tzuy*u.sY + w.X*usY_c2 + u.sTz + u.sX*u.sY2 + u.sX2*u.sY - u.sX2Y - u.sY2X;
				u.Tzux = w.Tzux + w.Y*u.sX + u.sY2;
				u.Tzuy = w.Tzuy + w.X*u.sY + u.sX2;
			}
			
			bool verify(Node &w){
				Node v = *this;
				testCalc(w);
				return memcmp(this, &v, sizeof(Node));
			}
		};
		
		void Xnormal(Node& u, Node& v, Node& o, Node& w){
			const score_t vX=v.X-o.X, vX2=(v.X*(v.X-1)-o.X*(o.X-1))/2, vY=v.Y, vY2=v.Y*(v.Y-1)/2, vZ=v.Z, vZ2=v.Z*(v.Z-1)/2, vQ=v.Q-o.Q, vTx=v.Tx, vTy=v.Ty-o.Ty, vTz=v.Tz-o.Tz;
			
			const score_t dusX_c2=((u.sX+vX)*(u.sX+vX-1)-u.sX*(u.sX-1))/2;
			u.X += vX;
			u.Q += vQ + w.Qux*vX + w.Quxy*vX*u.sY + w.Quxz*vX*u.sZ + w.Qux2*dusX_c2
				+ 2*w.Y*(vTy + vX*u.sZ2 + vX2*u.sZ - vZ2*vX - vX2*vZ) + 2*w.Z*(vTz + vX*u.sY2 + vX2*u.sY - vX2*vY - vY2*vX)
				+ 3*vX2*vY*vZ + 3*vY2*vX*vZ + 3*vZ2*vX*vY - 2*vX*u.sY2Z - 2*vX2*vZ*u.sY - vX2*u.sYZ - vX*vZ*u.sY2 - 2*vX*u.sZ2Y - 2*vX2*vY*u.sZ - vX*vY*u.sZ2 - 2*u.sY*vZ2*vX
				- 2*vY2*vX*u.sZ + 2*vX*u.sY*u.sZ2 + 2*vX*u.sY2*u.sZ + 2*vX2*u.sY*u.sZ
				+ 2*u.sTx*vX - 2*vTx*vX + 2*vTy*u.sY - 2*vTy*vY + 2*vTz*u.sZ - 2*vTz*vZ;
			u.Qux += w.Qux2*vX;
			u.Quy += 2*w.Z*vX2 + w.Quxy*vX + 2*vX*u.sZ2 + 2*vX2*u.sZ - 2*vX2*vZ - 2*vZ2*vX + 2*vTy;
			u.Quz += 2*w.Y*vX2 + w.Quxz*vX + 2*vX*u.sY2 + 2*vX2*u.sY - 2*vX2*vY - 2*vY2*vX + 2*vTz;
			u.Quyz += vX2;
			u.Quy2 += 2*w.Z*vX + 2*vX*u.sZ - vX*vZ;
			u.Quz2 += 2*w.Y*vX + 2*vX*u.sY - vX*vY;
			u.Ty += w.Tyux*vX + w.Z*dusX_c2 + vTy + vX*u.sZ2 + vX2*u.sZ - vZ2*vX - vX2*vZ;
			u.Tyux += w.Z*vX;
			u.Tyuz += vX2;
			u.Tz += w.Tzux*vX + w.Y*dusX_c2 + vTz + vX*u.sY2 + vX2*u.sY - vX2*vY - vY2*vX;
			u.Tzux += w.Y*vX;
			u.Tzuy += vX2;
			
			u.sX += vX;
			u.sXY += vX * vY;
			u.sXZ += vX * vZ;
			u.sX2 += vX2;
			u.sX2Y += vX2 * vY;
			u.sX2Z += vX2 * vZ;
			u.sY2X += vY2 * vX;
			u.sZ2X += vZ2 * vX;
			u.sX2YZ += vX2 * vY * vZ;
			u.sY2XZ += vY2 * vX * vZ;
			u.sZ2XY += vZ2 * vX * vY;
			u.sTxX += vTx * vX;
			u.sTy += vTy;
			u.sTyY += vTy * vY;
			u.sTz += vTz;
			u.sTzZ += vTz * vZ;
			u.sQ += vQ;
		}
		
		void Ynormal(Node& u, Node& v, Node& o, Node& w){
			const score_t vX=v.X, vX2=v.X*(v.X-1)/2, vY=v.Y-o.Y, vY2=(v.Y*(v.Y-1)-o.Y*(o.Y-1))/2, vZ=v.Z, vZ2=v.Z*(v.Z-1)/2, vQ=v.Q-o.Q, vTx=v.Tx-o.Tx, vTy=v.Ty, vTz=v.Tz-o.Tz;
			
			const score_t dusY_c2=((u.sY+vY)*(u.sY+vY-1)-u.sY*(u.sY-1))/2;
			u.Y += vY;
			u.Q += vQ + w.Quy*vY + w.Quxy*u.sX*vY + w.Quyz*vY*u.sZ + w.Quy2*dusY_c2
				+ 2*w.X*(vTx + vY*u.sZ2 + vY2*u.sZ - vZ2*vY - vY2*vZ) + 2*w.Z*(vTz + u.sX*vY2 + u.sX2*vY - vX2*vY - vY2*vX)
				+ 3*vX2*vY*vZ + 3*vY2*vX*vZ + 3*vZ2*vX*vY - 2*u.sX*vY2*vZ - 2*u.sX2Z*vY - u.sX2*vY*vZ - u.sXZ*vY2 - 2*u.sX*vZ2*vY - 2*vX2*vY*u.sZ - vX*vY*u.sZ2 - 2*vY*u.sZ2X
				- 2*vY2*vX*u.sZ + 2*u.sX*vY*u.sZ2 + 2*u.sX*vY2*u.sZ + 2*u.sX2*vY*u.sZ
				+ 2*vTx*u.sX - 2*vTx*vX + 2*u.sTy*vY - 2*vTy*vY + 2*vTz*u.sZ - 2*vTz*vZ;
			u.Qux += 2*w.Z*vY2 + w.Quxy*vY + 2*vY*u.sZ2 + 2*vY2*u.sZ - 2*vY2*vZ - 2*vZ2*vY + 2*vTx;
			u.Quy += w.Quy2*vY;
			u.Quz += 2*w.X*vY2 + w.Quyz*vY + 2*u.sX*vY2 + 2*u.sX2*vY - 2*vX2*vY - 2*vY2*vX + 2*vTz;
			u.Quxz += vY2;
			u.Qux2 += 2*w.Z*vY + 2*vY*u.sZ - vY*vZ;
			u.Quz2 += 2*w.X*vY + 2*u.sX*vY - vX*vY;
			u.Tx += w.Txuy*vY + w.Z*dusY_c2 + vTx + vY*u.sZ2 + vY2*u.sZ - vZ2*vY - vY2*vZ;
			u.Txuy += w.Z*vY;
			u.Txuz += vY2;
			u.Tz += w.Tzuy*vY + w.X*dusY_c2 + vTz + u.sX*vY2 + u.sX2*vY - vX2*vY - vY2*vX;
			u.Tzux += vY2;
			u.Tzuy += w.X*vY;
			
			u.sY += vY;
			u.sXY += vX * vY;
			u.sYZ += vY * vZ;
			u.sY2 += vY2;
			u.sX2Y += vX2 * vY;
			u.sY2X += vY2 * vX;
			u.sY2Z += vY2 * vZ;
			u.sZ2Y += vZ2 * vY;
			u.sX2YZ += vX2 * vY * vZ;
			u.sY2XZ += vY2 * vX * vZ;
			u.sZ2XY += vZ2 * vX * vY;
			u.sTx += vTx; 
			u.sTxX += vTx * vX;
			u.sTyY += vTy * vY;
			u.sTz += vTz;
			u.sTzZ += vTz * vZ;
			u.sQ += vQ;
		}
		
		void Znormal(Node& u, Node& v, Node& o, Node& w){
			const score_t vX=v.X, vX2=v.X*(v.X-1)/2, vY=v.Y, vY2=v.Y*(v.Y-1)/2, vZ=v.Z-o.Z, vZ2=(v.Z*(v.Z-1)-o.Z*(o.Z-1))/2, vQ=v.Q-o.Q, vTx=v.Tx-o.Tx, vTy=v.Ty-o.Ty, vTz=v.Tz;
			
			const score_t dusZ_c2=((u.sZ+vZ)*(u.sZ+vZ-1)-u.sZ*(u.sZ-1))/2;
			u.Z += vZ;
			u.Q += vQ + w.Quz*vZ + w.Quxz*u.sX*vZ + w.Quyz*u.sY*vZ + w.Quz2*dusZ_c2
				+ 2*w.X*(vTx + u.sY*vZ2 + u.sY2*vZ - vZ2*vY - vY2*vZ) + 2*w.Y*(vTy + u.sX*vZ2 + u.sX2*vZ - vZ2*vX - vX2*vZ)
				+ 3*vX2*vY*vZ + 3*vY2*vX*vZ + 3*vZ2*vX*vY - 2*u.sX*vY2*vZ - 2*vX2*vZ*u.sY - u.sX2*vY*vZ - vX*vZ*u.sY2 - 2*u.sX*vZ2*vY - 2*u.sX2Y*vZ - u.sXY*vZ2 - 2*u.sY*vZ2*vX
				- 2*u.sY2X*vZ + 2*u.sX*u.sY*vZ2 + 2*u.sX*u.sY2*vZ + 2*u.sX2*u.sY*vZ
				+ 2*vTx*u.sX - 2*vTx*vX + 2*vTy*u.sY - 2*vTy*vY + 2*u.sTz*vZ - 2*vTz*vZ;
			u.Qux += 2*w.Y*vZ2 + w.Quxz*vZ + 2*u.sY*vZ2 + 2*u.sY2*vZ - 2*vY2*vZ - 2*vZ2*vY + 2*vTx;
			u.Quy += 2*w.X*vZ2 + w.Quyz*vZ + 2*u.sX*vZ2 + 2*u.sX2*vZ - 2*vX2*vZ - 2*vZ2*vX + 2*vTy;
			u.Quz += w.Quz2*vZ;
			u.Quxy += vZ2;
			u.Qux2 += 2*w.Y*vZ + 2*u.sY*vZ - vY*vZ;
			u.Quy2 += 2*w.X*vZ + 2*u.sX*vZ - vX*vZ;
			u.Tx += w.Txuz*vZ + w.Y*dusZ_c2 + vTx + u.sY*vZ2 + u.sY2*vZ - vZ2*vY - vY2*vZ;
			u.Txuy += vZ2;
			u.Txuz += w.Y*vZ;
			u.Ty += w.Tyuz*vZ + w.X*dusZ_c2 + vTy + u.sX*vZ2 + u.sX2*vZ - vZ2*vX - vX2*vZ;
			u.Tyux += vZ2;
			u.Tyuz += w.X*vZ;
			
			u.sZ += vZ;
			u.sXZ += vX * vZ;
			u.sYZ += vY * vZ;
			u.sZ2 += vZ2;
			u.sX2Z += vX2 * vZ;
			u.sY2Z += vY2 * vZ;
			u.sZ2X += vZ2 * vX;
			u.sZ2Y += vZ2 * vY;
			u.sX2YZ += vX2 * vY * vZ;
			u.sY2XZ += vY2 * vX * vZ;
			u.sZ2XY += vZ2 * vX * vY;
			u.sTx += vTx; 
			u.sTxX += vTx * vX;
			u.sTy += vTy;
			u.sTyY += vTy * vY;
			u.sTzZ += vTz * vZ;
			u.sQ += vQ;
		}
		
		void Xspecial(Node& u, Node& w, Node& o){
			const score_t usY_c2 = u.sY*(u.sY-1)/2, usZ_c2 = u.sZ*(u.sZ-1)/2;
			const score_t wQ=w.Q-o.Q, wQux=w.Qux-o.Qux, wQuy=w.Quy-o.Quy, wQuz=w.Quz-o.Quz;
			const score_t wX=w.X-o.X, wQuy2=w.Quy2-o.Quy2, wQuz2=w.Quz2-o.Quz2, wQuyz=w.Quyz-o.Quyz;
			const score_t wTy=w.Ty-o.Ty, wTz=w.Tz-o.Tz, wTyux=w.Tyux-o.Tyux, wTyuz=w.Tyuz-o.Tyuz, wTzux=w.Tzux-o.Tzux, wTzuy=w.Tzuy-o.Tzuy;
			
			u.X += wX;
			u.Q += wQ + wQux*u.sX + wQuy*u.sY + wQuz*u.sZ + wQuyz*u.sY*u.sZ + wQuy2*usY_c2 + wQuz2*usZ_c2
				+ 2*wX*(u.sTx + u.sY*u.sZ2 + u.sY2*u.sZ - u.sZ2Y - u.sY2Z);
			u.Qux += wQux;
			u.Quy += 2*wX*u.sZ2 + wQuy + wQuy2*u.sY + wQuyz*u.sZ;
			u.Quz += 2*wX*u.sY2 + wQuz + wQuyz*u.sY + wQuz2*u.sZ;
			u.Quyz += wQuyz;
			u.Quy2 += 2*wX*u.sZ + wQuy2;
			u.Quz2 += 2*wX*u.sY + wQuz2;
			u.Ty += wTy + wTyux*u.sX + wTyuz*u.sZ + wX*usZ_c2;
			u.Tyux += wTyux;
			u.Tyuz += wTyuz + wX*u.sZ;
			u.Tz += wTz + wTzux*u.sX + wTzuy*u.sY + wX*usY_c2;
			u.Tzux += wTzux;
			u.Tzuy += wTzuy + wX*u.sY;
		}
		
		void Yspecial(Node& u, Node& w, Node& o){
			const score_t usX_c2 = u.sX*(u.sX-1)/2, usZ_c2 = u.sZ*(u.sZ-1)/2;
			const score_t wQ=w.Q-o.Q, wQux=w.Qux-o.Qux, wQuy=w.Quy-o.Quy, wQuz=w.Quz-o.Quz;
			const score_t wY=w.Y-o.Y, wQux2=w.Qux2-o.Qux2, wQuz2=w.Quz2-o.Quz2, wQuxz=w.Quxz-o.Quxz;
			const score_t wTx=w.Tx-o.Tx, wTz=w.Tz-o.Tz, wTxuy=w.Txuy-o.Txuy, wTxuz=w.Txuz-o.Txuz, wTzux=w.Tzux-o.Tzux, wTzuy=w.Tzuy-o.Tzuy;
			
			u.Y += wY;
			u.Q += wQ + wQux*u.sX + wQuy*u.sY + wQuz*u.sZ + wQuxz*u.sX*u.sZ + wQux2*usX_c2 + wQuz2*usZ_c2
				+ 2*wY*(u.sTy + u.sX*u.sZ2 + u.sX2*u.sZ - u.sZ2X - u.sX2Z);
			u.Qux += 2*wY*u.sZ2 + wQux + wQux2*u.sX + wQuxz*u.sZ;
			u.Quy += wQuy;
			u.Quz += 2*wY*u.sX2 + wQuz + wQuxz*u.sX + wQuz2*u.sZ;
			u.Quxz += wQuxz;
			u.Qux2 += 2*wY*u.sZ + wQux2;
			u.Quz2 += 2*wY*u.sX + wQuz2;
			u.Tx += wTx + wTxuy*u.sY + wTxuz*u.sZ + wY*usZ_c2;
			u.Txuy += wTxuy;
			u.Txuz += wTxuz + wY*u.sZ;
			u.Tz += wTz + wTzux*u.sX + wY*usX_c2 + wTzuy*u.sY;
			u.Tzux += wTzux + wY*u.sX;
			u.Tzuy += wTzuy;
		}
		
		void Zspecial(Node& u, Node& w, Node& o){
			const score_t usX_c2 = u.sX*(u.sX-1)/2, usY_c2 = u.sY*(u.sY-1)/2;
			const score_t wQ=w.Q-o.Q, wQux=w.Qux-o.Qux, wQuy=w.Quy-o.Quy, wQuz=w.Quz-o.Quz;
			const score_t wZ=w.Z-o.Z, wQux2=w.Qux2-o.Qux2, wQuy2=w.Quy2-o.Quy2, wQuxy=w.Quxy-o.Quxy;
			const score_t wTx=w.Tx-o.Tx, wTy=w.Ty-o.Ty, wTxuy=w.Txuy-o.Txuy, wTxuz=w.Txuz-o.Txuz, wTyux=w.Tyux-o.Tyux, wTyuz=w.Tyuz-o.Tyuz;
			
			u.Z += wZ;
			u.Q += wQ + wQux*u.sX + wQuy*u.sY + wQuz*u.sZ + wQuxy*u.sX*u.sY + wQux2*usX_c2 + wQuy2*usY_c2
				+ 2*wZ*(u.sTz + u.sX*u.sY2 + u.sX2*u.sY - u.sX2Y - u.sY2X);
			u.Qux += 2*wZ*u.sY2 + wQux + wQux2*u.sX + wQuxy*u.sY;
			u.Quy += 2*wZ*u.sX2 + wQuy + wQuxy*u.sX + wQuy2*u.sY;
			u.Quz += wQuz;
			u.Quxy += wQuxy;
			u.Qux2 += 2*wZ*u.sY + wQux2;
			u.Quy2 += 2*wZ*u.sX + wQuy2;
			u.Tx += wTx + wTxuy*u.sY + wZ*usY_c2 + wTxuz*u.sZ;
			u.Txuy += wTxuy + wZ*u.sY;
			u.Txuz += wTxuz;
			u.Ty += wTy + wTyux*u.sX + wZ*usX_c2 + wTyuz*u.sZ;
			u.Tyux += wTyux + wZ*u.sX;
			u.Tyuz += wTyuz; 
		}
		
		void doXnormal(Node& u, Node& v, Node& o, Node& w){
			const score_t vX=1, vX2=(v.X*(v.X-1)-o.X*(o.X-1))/2, vY=v.Y, vY2=v.Y*(v.Y-1)/2, vZ=v.Z, vZ2=v.Z*(v.Z-1)/2, vQ=v.Q-o.Q, vTx=v.Tx, vTy=v.Ty-o.Ty, vTz=v.Tz-o.Tz;
			
			const score_t dusX_c2=((u.sX+vX)*(u.sX+vX-1)-u.sX*(u.sX-1))/2;
			u.X += vX;
			u.Q += vQ + w.Qux*vX + w.Quxy*vX*u.sY + w.Quxz*vX*u.sZ + w.Qux2*dusX_c2
				+ 2*w.Y*(vTy + vX*u.sZ2 + vX2*u.sZ - vZ2*vX - vX2*vZ) + 2*w.Z*(vTz + vX*u.sY2 + vX2*u.sY - vX2*vY - vY2*vX)
				+ 3*vX2*vY*vZ + 3*vY2*vX*vZ + 3*vZ2*vX*vY - 2*vX*u.sY2Z - 2*vX2*vZ*u.sY - vX2*u.sYZ - vX*vZ*u.sY2 - 2*vX*u.sZ2Y - 2*vX2*vY*u.sZ - vX*vY*u.sZ2 - 2*u.sY*vZ2*vX
				- 2*vY2*vX*u.sZ + 2*vX*u.sY*u.sZ2 + 2*vX*u.sY2*u.sZ + 2*vX2*u.sY*u.sZ
				+ 2*u.sTx*vX - 2*vTx*vX + 2*vTy*u.sY - 2*vTy*vY + 2*vTz*u.sZ - 2*vTz*vZ;
			u.Qux += w.Qux2*vX;
			u.Quy += 2*w.Z*vX2 + w.Quxy*vX + 2*vX*u.sZ2 + 2*vX2*u.sZ - 2*vX2*vZ - 2*vZ2*vX + 2*vTy;
			u.Quz += 2*w.Y*vX2 + w.Quxz*vX + 2*vX*u.sY2 + 2*vX2*u.sY - 2*vX2*vY - 2*vY2*vX + 2*vTz;
			u.Quyz += vX2;
			u.Quy2 += 2*w.Z*vX + 2*vX*u.sZ - vX*vZ;
			u.Quz2 += 2*w.Y*vX + 2*vX*u.sY - vX*vY;
			u.Ty += w.Tyux*vX + w.Z*dusX_c2 + vTy + vX*u.sZ2 + vX2*u.sZ - vZ2*vX - vX2*vZ;
			u.Tyux += w.Z*vX;
			u.Tyuz += vX2;
			u.Tz += w.Tzux*vX + w.Y*dusX_c2 + vTz + vX*u.sY2 + vX2*u.sY - vX2*vY - vY2*vX;
			u.Tzux += w.Y*vX;
			u.Tzuy += vX2;
			
			u.sX += vX;
			u.sXY += vX * vY;
			u.sXZ += vX * vZ;
			u.sX2 += vX2;
			u.sX2Y += vX2 * vY;
			u.sX2Z += vX2 * vZ;
			u.sY2X += vY2 * vX;
			u.sZ2X += vZ2 * vX;
			u.sX2YZ += vX2 * vY * vZ;
			u.sY2XZ += vY2 * vX * vZ;
			u.sZ2XY += vZ2 * vX * vY;
			u.sTxX += vTx * vX;
			u.sTy += vTy;
			u.sTyY += vTy * vY;
			u.sTz += vTz;
			u.sTzZ += vTz * vZ;
			u.sQ += vQ;
		}
		
		void doYnormal(Node& u, Node& v, Node& o, Node& w){
			const score_t vX=v.X, vX2=v.X*(v.X-1)/2, vY=1, vY2=(v.Y*(v.Y-1)-o.Y*(o.Y-1))/2, vZ=v.Z, vZ2=v.Z*(v.Z-1)/2, vQ=v.Q-o.Q, vTx=v.Tx-o.Tx, vTy=v.Ty, vTz=v.Tz-o.Tz;
			
			const score_t dusY_c2=((u.sY+vY)*(u.sY+vY-1)-u.sY*(u.sY-1))/2;
			u.Y += vY;
			u.Q += vQ + w.Quy*vY + w.Quxy*u.sX*vY + w.Quyz*vY*u.sZ + w.Quy2*dusY_c2
				+ 2*w.X*(vTx + vY*u.sZ2 + vY2*u.sZ - vZ2*vY - vY2*vZ) + 2*w.Z*(vTz + u.sX*vY2 + u.sX2*vY - vX2*vY - vY2*vX)
				+ 3*vX2*vY*vZ + 3*vY2*vX*vZ + 3*vZ2*vX*vY - 2*u.sX*vY2*vZ - 2*u.sX2Z*vY - u.sX2*vY*vZ - u.sXZ*vY2 - 2*u.sX*vZ2*vY - 2*vX2*vY*u.sZ - vX*vY*u.sZ2 - 2*vY*u.sZ2X
				- 2*vY2*vX*u.sZ + 2*u.sX*vY*u.sZ2 + 2*u.sX*vY2*u.sZ + 2*u.sX2*vY*u.sZ
				+ 2*vTx*u.sX - 2*vTx*vX + 2*u.sTy*vY - 2*vTy*vY + 2*vTz*u.sZ - 2*vTz*vZ;
			u.Qux += 2*w.Z*vY2 + w.Quxy*vY + 2*vY*u.sZ2 + 2*vY2*u.sZ - 2*vY2*vZ - 2*vZ2*vY + 2*vTx;
			u.Quy += w.Quy2*vY;
			u.Quz += 2*w.X*vY2 + w.Quyz*vY + 2*u.sX*vY2 + 2*u.sX2*vY - 2*vX2*vY - 2*vY2*vX + 2*vTz;
			u.Quxz += vY2;
			u.Qux2 += 2*w.Z*vY + 2*vY*u.sZ - vY*vZ;
			u.Quz2 += 2*w.X*vY + 2*u.sX*vY - vX*vY;
			u.Tx += w.Txuy*vY + w.Z*dusY_c2 + vTx + vY*u.sZ2 + vY2*u.sZ - vZ2*vY - vY2*vZ;
			u.Txuy += w.Z*vY;
			u.Txuz += vY2;
			u.Tz += w.Tzuy*vY + w.X*dusY_c2 + vTz + u.sX*vY2 + u.sX2*vY - vX2*vY - vY2*vX;
			u.Tzux += vY2;
			u.Tzuy += w.X*vY;
			
			u.sY += vY;
			u.sXY += vX * vY;
			u.sYZ += vY * vZ;
			u.sY2 += vY2;
			u.sX2Y += vX2 * vY;
			u.sY2X += vY2 * vX;
			u.sY2Z += vY2 * vZ;
			u.sZ2Y += vZ2 * vY;
			u.sX2YZ += vX2 * vY * vZ;
			u.sY2XZ += vY2 * vX * vZ;
			u.sZ2XY += vZ2 * vX * vY;
			u.sTx += vTx; 
			u.sTxX += vTx * vX;
			u.sTyY += vTy * vY;
			u.sTz += vTz;
			u.sTzZ += vTz * vZ;
			u.sQ += vQ;
		}
		
		void doZnormal(Node& u, Node& v, Node& o, Node& w){
			const score_t vX=v.X, vX2=v.X*(v.X-1)/2, vY=v.Y, vY2=v.Y*(v.Y-1)/2, vZ=1, vZ2=(v.Z*(v.Z-1)-o.Z*(o.Z-1))/2, vQ=v.Q-o.Q, vTx=v.Tx-o.Tx, vTy=v.Ty-o.Ty, vTz=v.Tz;
			
			const score_t dusZ_c2=((u.sZ+vZ)*(u.sZ+vZ-1)-u.sZ*(u.sZ-1))/2;
			u.Z += vZ;
			u.Q += vQ + w.Quz*vZ + w.Quxz*u.sX*vZ + w.Quyz*u.sY*vZ + w.Quz2*dusZ_c2
				+ 2*w.X*(vTx + u.sY*vZ2 + u.sY2*vZ - vZ2*vY - vY2*vZ) + 2*w.Y*(vTy + u.sX*vZ2 + u.sX2*vZ - vZ2*vX - vX2*vZ)
				+ 3*vX2*vY*vZ + 3*vY2*vX*vZ + 3*vZ2*vX*vY - 2*u.sX*vY2*vZ - 2*vX2*vZ*u.sY - u.sX2*vY*vZ - vX*vZ*u.sY2 - 2*u.sX*vZ2*vY - 2*u.sX2Y*vZ - u.sXY*vZ2 - 2*u.sY*vZ2*vX
				- 2*u.sY2X*vZ + 2*u.sX*u.sY*vZ2 + 2*u.sX*u.sY2*vZ + 2*u.sX2*u.sY*vZ
				+ 2*vTx*u.sX - 2*vTx*vX + 2*vTy*u.sY - 2*vTy*vY + 2*u.sTz*vZ - 2*vTz*vZ;
			u.Qux += 2*w.Y*vZ2 + w.Quxz*vZ + 2*u.sY*vZ2 + 2*u.sY2*vZ - 2*vY2*vZ - 2*vZ2*vY + 2*vTx;
			u.Quy += 2*w.X*vZ2 + w.Quyz*vZ + 2*u.sX*vZ2 + 2*u.sX2*vZ - 2*vX2*vZ - 2*vZ2*vX + 2*vTy;
			u.Quz += w.Quz2*vZ;
			u.Quxy += vZ2;
			u.Qux2 += 2*w.Y*vZ + 2*u.sY*vZ - vY*vZ;
			u.Quy2 += 2*w.X*vZ + 2*u.sX*vZ - vX*vZ;
			u.Tx += w.Txuz*vZ + w.Y*dusZ_c2 + vTx + u.sY*vZ2 + u.sY2*vZ - vZ2*vY - vY2*vZ;
			u.Txuy += vZ2;
			u.Txuz += w.Y*vZ;
			u.Ty += w.Tyuz*vZ + w.X*dusZ_c2 + vTy + u.sX*vZ2 + u.sX2*vZ - vZ2*vX - vX2*vZ;
			u.Tyux += vZ2;
			u.Tyuz += w.X*vZ;
			
			u.sZ += vZ;
			u.sXZ += vX * vZ;
			u.sYZ += vY * vZ;
			u.sZ2 += vZ2;
			u.sX2Z += vX2 * vZ;
			u.sY2Z += vY2 * vZ;
			u.sZ2X += vZ2 * vX;
			u.sZ2Y += vZ2 * vY;
			u.sX2YZ += vX2 * vY * vZ;
			u.sY2XZ += vY2 * vX * vZ;
			u.sZ2XY += vZ2 * vX * vY;
			u.sTx += vTx; 
			u.sTxX += vTx * vX;
			u.sTy += vTy;
			u.sTyY += vTy * vY;
			u.sTzZ += vTz * vZ;
			u.sQ += vQ;
		}
		
		void doXspecial(Node& u, Node& w, Node& o){
			const score_t usY_c2 = u.sY*(u.sY-1)/2, usZ_c2 = u.sZ*(u.sZ-1)/2;
			const score_t wQ=w.Q-o.Q, wQux=w.Qux-o.Qux, wQuy=w.Quy-o.Quy, wQuz=w.Quz-o.Quz;
			const score_t wX=1, wQuy2=w.Quy2-o.Quy2, wQuz2=w.Quz2-o.Quz2, wQuyz=w.Quyz-o.Quyz;
			const score_t wTy=w.Ty-o.Ty, wTz=w.Tz-o.Tz, wTyux=w.Tyux-o.Tyux, wTyuz=w.Tyuz-o.Tyuz, wTzux=w.Tzux-o.Tzux, wTzuy=w.Tzuy-o.Tzuy;
			
			u.X += wX;
			u.Q += wQ + wQux*u.sX + wQuy*u.sY + wQuz*u.sZ + wQuyz*u.sY*u.sZ + wQuy2*usY_c2 + wQuz2*usZ_c2
				+ 2*wX*(u.sTx + u.sY*u.sZ2 + u.sY2*u.sZ - u.sZ2Y - u.sY2Z);
			u.Qux += wQux;
			u.Quy += 2*wX*u.sZ2 + wQuy + wQuy2*u.sY + wQuyz*u.sZ;
			u.Quz += 2*wX*u.sY2 + wQuz + wQuyz*u.sY + wQuz2*u.sZ;
			u.Quyz += wQuyz;
			u.Quy2 += 2*wX*u.sZ + wQuy2;
			u.Quz2 += 2*wX*u.sY + wQuz2;
			u.Ty += wTy + wTyux*u.sX + wTyuz*u.sZ + wX*usZ_c2;
			u.Tyux += wTyux;
			u.Tyuz += wTyuz + wX*u.sZ;
			u.Tz += wTz + wTzux*u.sX + wTzuy*u.sY + wX*usY_c2;
			u.Tzux += wTzux;
			u.Tzuy += wTzuy + wX*u.sY;
		}
		
		void doYspecial(Node& u, Node& w, Node& o){
			const score_t usX_c2 = u.sX*(u.sX-1)/2, usZ_c2 = u.sZ*(u.sZ-1)/2;
			const score_t wQ=w.Q-o.Q, wQux=w.Qux-o.Qux, wQuy=w.Quy-o.Quy, wQuz=w.Quz-o.Quz;
			const score_t wY=1, wQux2=w.Qux2-o.Qux2, wQuz2=w.Quz2-o.Quz2, wQuxz=w.Quxz-o.Quxz;
			const score_t wTx=w.Tx-o.Tx, wTz=w.Tz-o.Tz, wTxuy=w.Txuy-o.Txuy, wTxuz=w.Txuz-o.Txuz, wTzux=w.Tzux-o.Tzux, wTzuy=w.Tzuy-o.Tzuy;
			
			u.Y += wY;
			u.Q += wQ + wQux*u.sX + wQuy*u.sY + wQuz*u.sZ + wQuxz*u.sX*u.sZ + wQux2*usX_c2 + wQuz2*usZ_c2
				+ 2*wY*(u.sTy + u.sX*u.sZ2 + u.sX2*u.sZ - u.sZ2X - u.sX2Z);
			u.Qux += 2*wY*u.sZ2 + wQux + wQux2*u.sX + wQuxz*u.sZ;
			u.Quy += wQuy;
			u.Quz += 2*wY*u.sX2 + wQuz + wQuxz*u.sX + wQuz2*u.sZ;
			u.Quxz += wQuxz;
			u.Qux2 += 2*wY*u.sZ + wQux2;
			u.Quz2 += 2*wY*u.sX + wQuz2;
			u.Tx += wTx + wTxuy*u.sY + wTxuz*u.sZ + wY*usZ_c2;
			u.Txuy += wTxuy;
			u.Txuz += wTxuz + wY*u.sZ;
			u.Tz += wTz + wTzux*u.sX + wY*usX_c2 + wTzuy*u.sY;
			u.Tzux += wTzux + wY*u.sX;
			u.Tzuy += wTzuy;
		}
		
		void doZspecial(Node& u, Node& w, Node& o){
			const score_t usX_c2 = u.sX*(u.sX-1)/2, usY_c2 = u.sY*(u.sY-1)/2;
			const score_t wQ=w.Q-o.Q, wQux=w.Qux-o.Qux, wQuy=w.Quy-o.Quy, wQuz=w.Quz-o.Quz;
			const score_t wZ=1, wQux2=w.Qux2-o.Qux2, wQuy2=w.Quy2-o.Quy2, wQuxy=w.Quxy-o.Quxy;
			const score_t wTx=w.Tx-o.Tx, wTy=w.Ty-o.Ty, wTxuy=w.Txuy-o.Txuy, wTxuz=w.Txuz-o.Txuz, wTyux=w.Tyux-o.Tyux, wTyuz=w.Tyuz-o.Tyuz;
			
			u.Z += wZ;
			u.Q += wQ + wQux*u.sX + wQuy*u.sY + wQuz*u.sZ + wQuxy*u.sX*u.sY + wQux2*usX_c2 + wQuy2*usY_c2
				+ 2*wZ*(u.sTz + u.sX*u.sY2 + u.sX2*u.sY - u.sX2Y - u.sY2X);
			u.Qux += 2*wZ*u.sY2 + wQux + wQux2*u.sX + wQuxy*u.sY;
			u.Quy += 2*wZ*u.sX2 + wQuy + wQuxy*u.sX + wQuy2*u.sY;
			u.Quz += wQuz;
			u.Quxy += wQuxy;
			u.Qux2 += 2*wZ*u.sY + wQux2;
			u.Quy2 += 2*wZ*u.sX + wQuy2;
			u.Tx += wTx + wTxuy*u.sY + wZ*usY_c2 + wTxuz*u.sZ;
			u.Txuy += wTxuy + wZ*u.sY;
			u.Txuz += wTxuz;
			u.Ty += wTy + wTyux*u.sX + wZ*usX_c2 + wTyuz*u.sZ;
			u.Tyux += wTyux + wZ*u.sX;
			u.Tyuz += wTyuz; 
		}
		
		void undoXnormal(Node& u, Node& v, Node& o, Node& w){
			const score_t vX=-1, vX2=(v.X*(v.X-1)-o.X*(o.X-1))/2, vY=v.Y, vY2=v.Y*(v.Y-1)/2, vZ=v.Z, vZ2=v.Z*(v.Z-1)/2, vQ=v.Q-o.Q, vTx=v.Tx, vTy=v.Ty-o.Ty, vTz=v.Tz-o.Tz;
			
			const score_t dusX_c2=((u.sX+vX)*(u.sX+vX-1)-u.sX*(u.sX-1))/2;
			u.X += vX;
			u.Q += vQ + w.Qux*vX + w.Quxy*vX*u.sY + w.Quxz*vX*u.sZ + w.Qux2*dusX_c2
				+ 2*w.Y*(vTy + vX*u.sZ2 + vX2*u.sZ - vZ2*vX - vX2*vZ) + 2*w.Z*(vTz + vX*u.sY2 + vX2*u.sY - vX2*vY - vY2*vX)
				+ 3*vX2*vY*vZ + 3*vY2*vX*vZ + 3*vZ2*vX*vY - 2*vX*u.sY2Z - 2*vX2*vZ*u.sY - vX2*u.sYZ - vX*vZ*u.sY2 - 2*vX*u.sZ2Y - 2*vX2*vY*u.sZ - vX*vY*u.sZ2 - 2*u.sY*vZ2*vX
				- 2*vY2*vX*u.sZ + 2*vX*u.sY*u.sZ2 + 2*vX*u.sY2*u.sZ + 2*vX2*u.sY*u.sZ
				+ 2*u.sTx*vX - 2*vTx*vX + 2*vTy*u.sY - 2*vTy*vY + 2*vTz*u.sZ - 2*vTz*vZ;
			u.Qux += w.Qux2*vX;
			u.Quy += 2*w.Z*vX2 + w.Quxy*vX + 2*vX*u.sZ2 + 2*vX2*u.sZ - 2*vX2*vZ - 2*vZ2*vX + 2*vTy;
			u.Quz += 2*w.Y*vX2 + w.Quxz*vX + 2*vX*u.sY2 + 2*vX2*u.sY - 2*vX2*vY - 2*vY2*vX + 2*vTz;
			u.Quyz += vX2;
			u.Quy2 += 2*w.Z*vX + 2*vX*u.sZ - vX*vZ;
			u.Quz2 += 2*w.Y*vX + 2*vX*u.sY - vX*vY;
			u.Ty += w.Tyux*vX + w.Z*dusX_c2 + vTy + vX*u.sZ2 + vX2*u.sZ - vZ2*vX - vX2*vZ;
			u.Tyux += w.Z*vX;
			u.Tyuz += vX2;
			u.Tz += w.Tzux*vX + w.Y*dusX_c2 + vTz + vX*u.sY2 + vX2*u.sY - vX2*vY - vY2*vX;
			u.Tzux += w.Y*vX;
			u.Tzuy += vX2;
			
			u.sX += vX;
			u.sXY += vX * vY;
			u.sXZ += vX * vZ;
			u.sX2 += vX2;
			u.sX2Y += vX2 * vY;
			u.sX2Z += vX2 * vZ;
			u.sY2X += vY2 * vX;
			u.sZ2X += vZ2 * vX;
			u.sX2YZ += vX2 * vY * vZ;
			u.sY2XZ += vY2 * vX * vZ;
			u.sZ2XY += vZ2 * vX * vY;
			u.sTxX += vTx * vX;
			u.sTy += vTy;
			u.sTyY += vTy * vY;
			u.sTz += vTz;
			u.sTzZ += vTz * vZ;
			u.sQ += vQ;
		}
		
		void undoYnormal(Node& u, Node& v, Node& o, Node& w){
			const score_t vX=v.X, vX2=v.X*(v.X-1)/2, vY=-1, vY2=(v.Y*(v.Y-1)-o.Y*(o.Y-1))/2, vZ=v.Z, vZ2=v.Z*(v.Z-1)/2, vQ=v.Q-o.Q, vTx=v.Tx-o.Tx, vTy=v.Ty, vTz=v.Tz-o.Tz;
			
			const score_t dusY_c2=((u.sY+vY)*(u.sY+vY-1)-u.sY*(u.sY-1))/2;
			u.Y += vY;
			u.Q += vQ + w.Quy*vY + w.Quxy*u.sX*vY + w.Quyz*vY*u.sZ + w.Quy2*dusY_c2
				+ 2*w.X*(vTx + vY*u.sZ2 + vY2*u.sZ - vZ2*vY - vY2*vZ) + 2*w.Z*(vTz + u.sX*vY2 + u.sX2*vY - vX2*vY - vY2*vX)
				+ 3*vX2*vY*vZ + 3*vY2*vX*vZ + 3*vZ2*vX*vY - 2*u.sX*vY2*vZ - 2*u.sX2Z*vY - u.sX2*vY*vZ - u.sXZ*vY2 - 2*u.sX*vZ2*vY - 2*vX2*vY*u.sZ - vX*vY*u.sZ2 - 2*vY*u.sZ2X
				- 2*vY2*vX*u.sZ + 2*u.sX*vY*u.sZ2 + 2*u.sX*vY2*u.sZ + 2*u.sX2*vY*u.sZ
				+ 2*vTx*u.sX - 2*vTx*vX + 2*u.sTy*vY - 2*vTy*vY + 2*vTz*u.sZ - 2*vTz*vZ;
			u.Qux += 2*w.Z*vY2 + w.Quxy*vY + 2*vY*u.sZ2 + 2*vY2*u.sZ - 2*vY2*vZ - 2*vZ2*vY + 2*vTx;
			u.Quy += w.Quy2*vY;
			u.Quz += 2*w.X*vY2 + w.Quyz*vY + 2*u.sX*vY2 + 2*u.sX2*vY - 2*vX2*vY - 2*vY2*vX + 2*vTz;
			u.Quxz += vY2;
			u.Qux2 += 2*w.Z*vY + 2*vY*u.sZ - vY*vZ;
			u.Quz2 += 2*w.X*vY + 2*u.sX*vY - vX*vY;
			u.Tx += w.Txuy*vY + w.Z*dusY_c2 + vTx + vY*u.sZ2 + vY2*u.sZ - vZ2*vY - vY2*vZ;
			u.Txuy += w.Z*vY;
			u.Txuz += vY2;
			u.Tz += w.Tzuy*vY + w.X*dusY_c2 + vTz + u.sX*vY2 + u.sX2*vY - vX2*vY - vY2*vX;
			u.Tzux += vY2;
			u.Tzuy += w.X*vY;
			
			u.sY += vY;
			u.sXY += vX * vY;
			u.sYZ += vY * vZ;
			u.sY2 += vY2;
			u.sX2Y += vX2 * vY;
			u.sY2X += vY2 * vX;
			u.sY2Z += vY2 * vZ;
			u.sZ2Y += vZ2 * vY;
			u.sX2YZ += vX2 * vY * vZ;
			u.sY2XZ += vY2 * vX * vZ;
			u.sZ2XY += vZ2 * vX * vY;
			u.sTx += vTx; 
			u.sTxX += vTx * vX;
			u.sTyY += vTy * vY;
			u.sTz += vTz;
			u.sTzZ += vTz * vZ;
			u.sQ += vQ;
		}
		
		void undoZnormal(Node& u, Node& v, Node& o, Node& w){
			const score_t vX=v.X, vX2=v.X*(v.X-1)/2, vY=v.Y, vY2=v.Y*(v.Y-1)/2, vZ=-1, vZ2=(v.Z*(v.Z-1)-o.Z*(o.Z-1))/2, vQ=v.Q-o.Q, vTx=v.Tx-o.Tx, vTy=v.Ty-o.Ty, vTz=v.Tz;
			
			const score_t dusZ_c2=((u.sZ+vZ)*(u.sZ+vZ-1)-u.sZ*(u.sZ-1))/2;
			u.Z += vZ;
			u.Q += vQ + w.Quz*vZ + w.Quxz*u.sX*vZ + w.Quyz*u.sY*vZ + w.Quz2*dusZ_c2
				+ 2*w.X*(vTx + u.sY*vZ2 + u.sY2*vZ - vZ2*vY - vY2*vZ) + 2*w.Y*(vTy + u.sX*vZ2 + u.sX2*vZ - vZ2*vX - vX2*vZ)
				+ 3*vX2*vY*vZ + 3*vY2*vX*vZ + 3*vZ2*vX*vY - 2*u.sX*vY2*vZ - 2*vX2*vZ*u.sY - u.sX2*vY*vZ - vX*vZ*u.sY2 - 2*u.sX*vZ2*vY - 2*u.sX2Y*vZ - u.sXY*vZ2 - 2*u.sY*vZ2*vX
				- 2*u.sY2X*vZ + 2*u.sX*u.sY*vZ2 + 2*u.sX*u.sY2*vZ + 2*u.sX2*u.sY*vZ
				+ 2*vTx*u.sX - 2*vTx*vX + 2*vTy*u.sY - 2*vTy*vY + 2*u.sTz*vZ - 2*vTz*vZ;
			u.Qux += 2*w.Y*vZ2 + w.Quxz*vZ + 2*u.sY*vZ2 + 2*u.sY2*vZ - 2*vY2*vZ - 2*vZ2*vY + 2*vTx;
			u.Quy += 2*w.X*vZ2 + w.Quyz*vZ + 2*u.sX*vZ2 + 2*u.sX2*vZ - 2*vX2*vZ - 2*vZ2*vX + 2*vTy;
			u.Quz += w.Quz2*vZ;
			u.Quxy += vZ2;
			u.Qux2 += 2*w.Y*vZ + 2*u.sY*vZ - vY*vZ;
			u.Quy2 += 2*w.X*vZ + 2*u.sX*vZ - vX*vZ;
			u.Tx += w.Txuz*vZ + w.Y*dusZ_c2 + vTx + u.sY*vZ2 + u.sY2*vZ - vZ2*vY - vY2*vZ;
			u.Txuy += vZ2;
			u.Txuz += w.Y*vZ;
			u.Ty += w.Tyuz*vZ + w.X*dusZ_c2 + vTy + u.sX*vZ2 + u.sX2*vZ - vZ2*vX - vX2*vZ;
			u.Tyux += vZ2;
			u.Tyuz += w.X*vZ;
			
			u.sZ += vZ;
			u.sXZ += vX * vZ;
			u.sYZ += vY * vZ;
			u.sZ2 += vZ2;
			u.sX2Z += vX2 * vZ;
			u.sY2Z += vY2 * vZ;
			u.sZ2X += vZ2 * vX;
			u.sZ2Y += vZ2 * vY;
			u.sX2YZ += vX2 * vY * vZ;
			u.sY2XZ += vY2 * vX * vZ;
			u.sZ2XY += vZ2 * vX * vY;
			u.sTx += vTx; 
			u.sTxX += vTx * vX;
			u.sTy += vTy;
			u.sTyY += vTy * vY;
			u.sTzZ += vTz * vZ;
			u.sQ += vQ;
		}
		
		void undoXspecial(Node& u, Node& w, Node& o){
			const score_t usY_c2 = u.sY*(u.sY-1)/2, usZ_c2 = u.sZ*(u.sZ-1)/2;
			const score_t wQ=w.Q-o.Q, wQux=w.Qux-o.Qux, wQuy=w.Quy-o.Quy, wQuz=w.Quz-o.Quz;
			const score_t wX=-1, wQuy2=w.Quy2-o.Quy2, wQuz2=w.Quz2-o.Quz2, wQuyz=w.Quyz-o.Quyz;
			const score_t wTy=w.Ty-o.Ty, wTz=w.Tz-o.Tz, wTyux=w.Tyux-o.Tyux, wTyuz=w.Tyuz-o.Tyuz, wTzux=w.Tzux-o.Tzux, wTzuy=w.Tzuy-o.Tzuy;
			
			u.X += wX;
			u.Q += wQ + wQux*u.sX + wQuy*u.sY + wQuz*u.sZ + wQuyz*u.sY*u.sZ + wQuy2*usY_c2 + wQuz2*usZ_c2
				+ 2*wX*(u.sTx + u.sY*u.sZ2 + u.sY2*u.sZ - u.sZ2Y - u.sY2Z);
			u.Qux += wQux;
			u.Quy += 2*wX*u.sZ2 + wQuy + wQuy2*u.sY + wQuyz*u.sZ;
			u.Quz += 2*wX*u.sY2 + wQuz + wQuyz*u.sY + wQuz2*u.sZ;
			u.Quyz += wQuyz;
			u.Quy2 += 2*wX*u.sZ + wQuy2;
			u.Quz2 += 2*wX*u.sY + wQuz2;
			u.Ty += wTy + wTyux*u.sX + wTyuz*u.sZ + wX*usZ_c2;
			u.Tyux += wTyux;
			u.Tyuz += wTyuz + wX*u.sZ;
			u.Tz += wTz + wTzux*u.sX + wTzuy*u.sY + wX*usY_c2;
			u.Tzux += wTzux;
			u.Tzuy += wTzuy + wX*u.sY;
		}
		
		void undoYspecial(Node& u, Node& w, Node& o){
			const score_t usX_c2 = u.sX*(u.sX-1)/2, usZ_c2 = u.sZ*(u.sZ-1)/2;
			const score_t wQ=w.Q-o.Q, wQux=w.Qux-o.Qux, wQuy=w.Quy-o.Quy, wQuz=w.Quz-o.Quz;
			const score_t wY=-1, wQux2=w.Qux2-o.Qux2, wQuz2=w.Quz2-o.Quz2, wQuxz=w.Quxz-o.Quxz;
			const score_t wTx=w.Tx-o.Tx, wTz=w.Tz-o.Tz, wTxuy=w.Txuy-o.Txuy, wTxuz=w.Txuz-o.Txuz, wTzux=w.Tzux-o.Tzux, wTzuy=w.Tzuy-o.Tzuy;
			
			u.Y += wY;
			u.Q += wQ + wQux*u.sX + wQuy*u.sY + wQuz*u.sZ + wQuxz*u.sX*u.sZ + wQux2*usX_c2 + wQuz2*usZ_c2
				+ 2*wY*(u.sTy + u.sX*u.sZ2 + u.sX2*u.sZ - u.sZ2X - u.sX2Z);
			u.Qux += 2*wY*u.sZ2 + wQux + wQux2*u.sX + wQuxz*u.sZ;
			u.Quy += wQuy;
			u.Quz += 2*wY*u.sX2 + wQuz + wQuxz*u.sX + wQuz2*u.sZ;
			u.Quxz += wQuxz;
			u.Qux2 += 2*wY*u.sZ + wQux2;
			u.Quz2 += 2*wY*u.sX + wQuz2;
			u.Tx += wTx + wTxuy*u.sY + wTxuz*u.sZ + wY*usZ_c2;
			u.Txuy += wTxuy;
			u.Txuz += wTxuz + wY*u.sZ;
			u.Tz += wTz + wTzux*u.sX + wY*usX_c2 + wTzuy*u.sY;
			u.Tzux += wTzux + wY*u.sX;
			u.Tzuy += wTzuy;
		}
		
		void undoZspecial(Node& u, Node& w, Node& o){
			const score_t usX_c2 = u.sX*(u.sX-1)/2, usY_c2 = u.sY*(u.sY-1)/2;
			const score_t wQ=w.Q-o.Q, wQux=w.Qux-o.Qux, wQuy=w.Quy-o.Quy, wQuz=w.Quz-o.Quz;
			const score_t wZ=-1, wQux2=w.Qux2-o.Qux2, wQuy2=w.Quy2-o.Quy2, wQuxy=w.Quxy-o.Quxy;
			const score_t wTx=w.Tx-o.Tx, wTy=w.Ty-o.Ty, wTxuy=w.Txuy-o.Txuy, wTxuz=w.Txuz-o.Txuz, wTyux=w.Tyux-o.Tyux, wTyuz=w.Tyuz-o.Tyuz;
			
			u.Z += wZ;
			u.Q += wQ + wQux*u.sX + wQuy*u.sY + wQuz*u.sZ + wQuxy*u.sX*u.sY + wQux2*usX_c2 + wQuy2*usY_c2
				+ 2*wZ*(u.sTz + u.sX*u.sY2 + u.sX2*u.sY - u.sX2Y - u.sY2X);
			u.Qux += 2*wZ*u.sY2 + wQux + wQux2*u.sX + wQuxy*u.sY;
			u.Quy += 2*wZ*u.sX2 + wQuy + wQuxy*u.sX + wQuy2*u.sY;
			u.Quz += wQuz;
			u.Quxy += wQuxy;
			u.Qux2 += 2*wZ*u.sY + wQux2;
			u.Quy2 += 2*wZ*u.sX + wQuy2;
			u.Tx += wTx + wTxuy*u.sY + wZ*usY_c2 + wTxuz*u.sZ;
			u.Txuy += wTxuy + wZ*u.sY;
			u.Txuz += wTxuz;
			u.Ty += wTy + wTyux*u.sX + wZ*usX_c2 + wTyuz*u.sZ;
			u.Tyux += wTyux + wZ*u.sX;
			u.Tyuz += wTyuz; 
		}
		
		
		vector<vector<int> > leafParent;
		score_t totalScore = 0;
		Node dummy, nodeX, nodeY, nodeZ;
		vector<vector<Node> > nodes, nodesTotal;
		vector<int> roots;
		vector<int> color;
		
		int buildTree(int v, vector<TripartitionInitializer::Node> &iNodes, vector<vector<int> > &children){
			for (int c: children[v]){
				int u = buildTree(c, iNodes, children);
				iNodes[u].up = v;
				if (iNodes[u].maxHeight + 1 > iNodes[u].maxHeight){
					int t = iNodes[v].maxChild;
					iNodes[v].maxChild = u;
					iNodes[v].maxHeight = iNodes[u].maxHeight + 1;
					u = t;
				}
				if (u != -1 && iNodes[u].max2Height + 1 > iNodes[u].maxHeight){
					iNodes[v].max2Child = u;
					iNodes[v].max2Height = iNodes[u].maxHeight + 1;
				}
			}
			int u = iNodes[v].maxChild;
			if (u != -1 && iNodes[u].maxHeight >= iNodes[u].path + 1 + iNodes[v].max2Height){
				iNodes[u].up = -1;
				if (iNodes[u].downMost == -1){
					iNodes[v].up = u;
					iNodes[u].down = v;
				}
				else {
					iNodes[v].up = iNodes[u].downMost;
					iNodes[iNodes[u].downMost].down = v;
				}
				iNodes[u].downMost = v;
				iNodes[u].path += 1 + iNodes[v].max2Height;
				return u;
			}
			else return v;
		}
		
		Partition(const TripartitionInitializer &init, int k): nodes(init.roots[k].size()), leafParent(init.roots[k].size(), vector<int>(init.nTaxa, -1)),
			nodesTotal(init.roots[k].size()), color(init.nTaxa, -1){
			nodeX.X = 1;
			nodeY.Y = 1;
			nodeZ.Z = 1;
			
			//cerr << "here1!" << endl;
			for (int i = 0; i < init.roots[k].size(); i++){
				int nInternal = 0;
				for (const pair<int, int> pc: init.parentChild[k][i]){
					if (pc.first - init.nTaxa + 1 > nInternal) nInternal = pc.first - init.nTaxa + 1;
				}
				vector<TripartitionInitializer::Node> iNodes(nInternal);
				vector<vector<int> > children(nInternal);
				for (const pair<int, int> pc: init.parentChild[k][i]){
					if (pc.second < init.nTaxa) leafParent[i][pc.second] = pc.first - init.nTaxa;
					else children[pc.first - init.nTaxa].push_back(pc.second - init.nTaxa);
				}
				//cerr << init.roots[i] - init.nTaxa << endl;
				/*for (int i = 0; i < children.size(); i++){
					for (int j: children[i]){
						cerr << i << " " << j << endl;
					}
				}*/
				//for (int j = 0; j < init.nTaxa; j++) cerr << j << "->" << leafParent[i][j] << endl;
				roots.push_back(buildTree(init.roots[k][i] - init.nTaxa, iNodes, children));
				for (int j = 0; j < nInternal; j++){
					nodes[i].emplace_back();
					nodes[i][j].up = iNodes[j].up;
					nodes[i][j].down = iNodes[j].down;
					nodesTotal[i].emplace_back();
					nodesTotal[i][j].up = iNodes[j].up;
					nodesTotal[i][j].down = iNodes[j].down;
					//cerr << j << " " << nodes[i][j].up << " " << nodes[i][j].down << endl;
				}
			}
			//cerr << "Init finished\n";
		}
		
		void add(int x, int i){
			//cerr << nodes[0][0].X << " " << nodes[0][0].Y << " " << nodes[0][0].Z << ":" << nodes[0][0].Q << "|" << nodes[0][1].X << " " << nodes[0][1].Y << " " << nodes[0][1].Z << ":" << nodes[0][1].Q << "|" << nodes[0][2].X << " " << nodes[0][2].Y << " " << nodes[0][2].Z << ":" << nodes[0][2].Q << "|" << nodes[0][3].X << " " << nodes[0][3].Y << " " << nodes[0][3].Z << ":" << nodes[0][3].Q << endl << endl;
			totalScore = 0;
			for (int k = 0; k < leafParent.size(); k++){
				int u = leafParent[k][i], v = -2;
				Node *pv = (x == 0) ? &nodeZ : (x == 1) ? &nodeX : &nodeY;
				Node vOld = dummy;
				if (u == -1){
					totalScore += nodes[k][roots[k]].Q;
				}
				while (u != -1){
					Node *pu = &nodes[k][u];
					Node uOld = *pu;
					if (pu->down == v){
						if (x == 0) doZspecial(*pu, *pv, vOld);//{undoZspecial(*pu, vOld); doZspecial(*pu, *pv);}
						else if (x == 1) doXspecial(*pu, *pv, vOld);//{undoXspecial(*pu, vOld); doXspecial(*pu, *pv);}
						else doYspecial(*pu, *pv, vOld);//{undoYspecial(*pu, vOld); doYspecial(*pu, *pv);}
						
						//if (pu->verify(*pv)) cerr << "ERROR!";
					}
					else {
						Node &w = (pu->down == -1) ? dummy : nodes[k][pu->down];
						if (x == 0) doZnormal(*pu, *pv, vOld, w);//{undoZnormal(*pu, vOld, w); doZnormal(*pu, *pv, w);}
						else if (x == 1) doXnormal(*pu, *pv, vOld, w);//{undoXnormal(*pu, vOld, w); doXnormal(*pu, *pv, w);}
						else doYnormal(*pu, *pv, vOld, w);//{undoYnormal(*pu, vOld, w); doYnormal(*pu, *pv, w);}
						
						//if (pu->verify(w)) cerr << "ERROR!";
					}
					vOld = uOld;
					v = u;
					pv = pu;
					u = pu->up;
				}
				totalScore += pv->Q;
			}
			//cerr << nodes[0][0].X << " " << nodes[0][0].Y << " " << nodes[0][0].Z << ":" << nodes[0][0].Q << "|" << nodes[0][1].X << " " << nodes[0][1].Y << " " << nodes[0][1].Z << ":" << nodes[0][1].Q << "|" << nodes[0][2].X << " " << nodes[0][2].Y << " " << nodes[0][2].Z << ":" << nodes[0][2].Q << "|" << nodes[0][3].X << " " << nodes[0][3].Y << " " << nodes[0][3].Z << ":" << nodes[0][3].Q << endl << endl;
		}
		
		void rmv(int x, int i){
			totalScore = 0;
			for (int k = 0; k < leafParent.size(); k++){
				int u = leafParent[k][i], v = -2;
				Node *pv = &dummy;
				Node vOld = (x == 0) ? nodeZ : (x == 1) ? nodeX : nodeY;
				if (u == -1){
					totalScore += nodes[k][roots[k]].Q;
				}
				while (u != -1){
					Node *pu = &nodes[k][u];
					Node uOld = *pu;
					if (pu->down == v){
						if (x == 0) undoZspecial(*pu, *pv, vOld);//{undoZspecial(*pu, vOld); doZspecial(*pu, *pv);}
						else if (x == 1) undoXspecial(*pu, *pv, vOld);//{undoXspecial(*pu, vOld); doXspecial(*pu, *pv);}
						else undoYspecial(*pu, *pv, vOld);//{undoYspecial(*pu, vOld); doYspecial(*pu, *pv);}
					}
					else {
						Node &w = (pu->down == -1) ? dummy : nodes[k][pu->down];
						if (x == 0) undoZnormal(*pu, *pv, vOld, w);//{undoZnormal(*pu, vOld, w); doZnormal(*pu, *pv, w);}
						else if (x == 1) undoXnormal(*pu, *pv, vOld, w);//{undoXnormal(*pu, vOld, w); doXnormal(*pu, *pv, w);}
						else undoYnormal(*pu, *pv, vOld, w);//{undoYnormal(*pu, vOld, w); doYnormal(*pu, *pv, w);}
					}
					vOld = uOld;
					v = u;
					pv = pu;
					u = pu->up;
				}
				totalScore += pv->Q;
			}
		}
		
		void update(int x, int i){
			int y = color[i];
			if (x == y) return;
			color[i] = x;
			if (y == -1){
				add(x, i);
				return;
			}
			if (x == -1){
				rmv(y, i);
				return;
			}
			totalScore = 0;
			for (int k = 0; k < leafParent.size(); k++){
				int u = leafParent[k][i], v = -2;
				Node *pv = (x == 0) ? &nodeZ : (x == 1) ? &nodeX : &nodeY;
				Node vMid = dummy;
				Node vOld = (y == 0) ? nodeZ : (y == 1) ? nodeX : nodeY;
				if (u == -1){
					totalScore += nodes[k][roots[k]].Q;
				}
				while (u != -1){
					Node *pu = &nodes[k][u];
					Node uOld = *pu, uMid;
					if (pu->down == v){
						if (y == 0) undoZspecial(*pu, vMid, vOld);
						else if (y == 1) undoXspecial(*pu, vMid, vOld);
						else undoYspecial(*pu, vMid, vOld);
						uMid = *pu;
						if (x == 0) doZspecial(*pu, *pv, vMid);
						else if (x == 1) doXspecial(*pu, *pv, vMid);
						else doYspecial(*pu, *pv, vMid);
					}
					else {
						Node &w = (pu->down == -1) ? dummy : nodes[k][pu->down];
						if (y == 0) undoZnormal(*pu, vMid, vOld, w);
						else if (y == 1) undoXnormal(*pu, vMid, vOld, w);
						else undoYnormal(*pu, vMid, vOld, w);
						uMid = *pu;
						if (x == 0) doZnormal(*pu, *pv, vMid, w);
						else if (x == 1) doXnormal(*pu, *pv, vMid, w);
						else doYnormal(*pu, *pv, vMid, w);
					}
					vOld = uOld;
					vMid = uMid;
					v = u;
					pv = pu;
					u = pu->up;
				}
				totalScore += pv->Q;
			}
		}
		
		score_t score(){
			return totalScore;
		}
	};
	
	vector<Partition> parts;
	
	Tripartition(const TripartitionInitializer &init){
		for (int p = 0; p < init.roots.size(); p++) parts.emplace_back(init, p);
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

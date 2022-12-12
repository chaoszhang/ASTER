#define OBJECTIVE_VERSION "0"

#include<vector>
#include<array>
#include<thread>
#include<future>
#include<bitset>
#include<memory>
#include "threadpool.hpp"

//#define G_SUPPORT

using namespace std;

inline long long XXYY(long long x0, long long x1, long long x2, long long y0, long long y1, long long y2){
	return x0 * (x0 - 1) * y1 * y2 + x1 * (x1 - 1) * y2 * y0 + x2 * (x2 - 1) * y0 * y1
	     + y0 * (y0 - 1) * x1 * x2 + y1 * (y1 - 1) * x2 * x0 + y2 * (y2 - 1) * x0 * x1;
}

inline score_t scorePos(const array<array<unsigned short, 4>, 3> &cnt, const array<float, 4> &pi){
//lst = simplify([sABCD(R, R, Y, Y); sABCD(A, A, Y, Y); sABCD(C, C, R, R); sABCD(A, A, C, C)])
//sol = [Pa^2*Pc^2; -Pc^2*(Pa + Pg)^2; -Pa^2*(Pc + Pt)^2; (Pa + Pg)^2*(Pc + Pt)^2]
//lst2 = simplify([sABCD(Y, Y, R, R); sABCD(Y, Y, A, A); sABCD(R, R, C, C); sABCD(C, C, A, A)])

// (Pa^2+Pg^2)*(Pc^2+Pt^2)*sABCD(R, R, Y, Y)
// -Pr^2*(Pc^2+Pt^2)*sABCD(A, A, Y, Y) -Pr^2*(Pc^2+Pt^2)*sABCD(G, G, Y, Y) -(Pa^2+Pg^2)*Py^2*sABCD(C, C, R, R) -(Pa^2+Pg^2)*Py^2*sABCD(T, T, R, R)
// +Pr^2*Py^2*sABCD(A, A, C, C) +Pr^2*Py^2*sABCD(G, G, C, C) +Pr^2*Py^2*sABCD(A, A, T, T) +Pr^2*Py^2*sABCD(G, G, T, T)

	const score_t A = pi[0], C = pi[1], G = pi[2], T = pi[3];
	const score_t R = A + G, Y = C + T, R2 = A * A + G * G, Y2 = C * C + T * T;
	const long long a0 = cnt[0][0], c0 = cnt[0][1], g0 = cnt[0][2], t0 = cnt[0][3], r0 = a0 + g0, y0 = c0 + t0;
	const long long a1 = cnt[1][0], c1 = cnt[1][1], g1 = cnt[1][2], t1 = cnt[1][3], r1 = a1 + g1, y1 = c1 + t1;
	const long long a2 = cnt[2][0], c2 = cnt[2][1], g2 = cnt[2][2], t2 = cnt[2][3], r2 = a2 + g2, y2 = c2 + t2;
	
	const long long rryy = XXYY(r0, r1, r2, y0, y1, y2);

	const long long aayy = XXYY(a0, a1, a2, y0, y1, y2);
	const long long ggyy = XXYY(g0, g1, g2, y0, y1, y2);
	const long long rrcc = XXYY(r0, r1, r2, c0, c1, c2);
	const long long rrtt = XXYY(r0, r1, r2, t0, t1, t2);
	
	const long long aacc = XXYY(a0, a1, a2, c0, c1, c2);
	const long long aatt = XXYY(a0, a1, a2, t0, t1, t2);
	const long long ggcc = XXYY(g0, g1, g2, c0, c1, c2);
	const long long ggtt = XXYY(g0, g1, g2, t0, t1, t2);
	
	return rryy * R2 * Y2 - (aayy + ggyy) * (R * R) * Y2 - (rrcc + rrtt) * R2 * (Y * Y)
	     + (aacc + aatt + ggcc + ggtt) * (R * R) * (Y * Y);
}

const int SEQ_BATCH_BIT = 25, SEQ_BATCH_SIZE = 1 << SEQ_BATCH_BIT;

struct TripartitionInitializer{
	struct Sequence{
		vector<shared_ptr<bitset<SEQ_BATCH_SIZE> > > seq;
		size_t p = 0;

		Sequence(){
			seq.emplace_back(new bitset<SEQ_BATCH_SIZE>());
		}

		void append(char c) {
			bitset<SEQ_BATCH_SIZE>& b = *(seq.back());
			switch (c) {
				case 'A': case 'a': b.set(p++, 1); b.set(p++, 0); b.set(p++, 0); b.set(p++, 0); break;
				case 'C': case 'c': b.set(p++, 0); b.set(p++, 1); b.set(p++, 0); b.set(p++, 0); break;
				case 'G': case 'g': b.set(p++, 0); b.set(p++, 0); b.set(p++, 1); b.set(p++, 0); break;
				case 'T': case 't': b.set(p++, 0); b.set(p++, 0); b.set(p++, 0); b.set(p++, 1); break;
				default: b.set(p++, 0); b.set(p++, 0); b.set(p++, 0); b.set(p++, 0);
			}
			if (p == SEQ_BATCH_SIZE) {
				p = 0;
				seq.emplace_back(new bitset<SEQ_BATCH_SIZE>());
			}
		}

		size_t len() const {
			return ((((size_t) seq.size() - 1) << SEQ_BATCH_BIT) + p) >> 2;
		}

		bool get(size_t id, size_t pos) const {
			size_t i = id * 4 + pos;
			bitset<SEQ_BATCH_SIZE>& b = *(seq[i >> SEQ_BATCH_BIT]);
			return b[i & (SEQ_BATCH_SIZE - 1)];
		}

		char get(size_t id) const {
			if (get(id, 0)) return 'A';
			if (get(id, 1)) return 'C';
			if (get(id, 2)) return 'G';
			if (get(id, 3)) return 'T';
			return 'N';
		}

		void add(size_t id, array<size_t, 4>& cnt) const {
			size_t i = id * 4;
			bitset<SEQ_BATCH_SIZE>& b = *(seq[i >> SEQ_BATCH_BIT]);
			i &= SEQ_BATCH_SIZE - 1;
			for (int t = 0; t < 4; t++) {
				cnt[t] += b[i++];
			}
		}

		inline void add(size_t id, array<unsigned short, 4>& cnt) const {
			size_t i = id * 4;
			bitset<SEQ_BATCH_SIZE>& b = *(seq[i >> SEQ_BATCH_BIT]);
			i &= SEQ_BATCH_SIZE - 1;
			for (int t = 0; t < 4; t++) {
				cnt[t] += b[i++];
			}
		}

		inline void rmv(size_t id, array<unsigned short, 4>& cnt) const {
			size_t i = id * 4;
			bitset<SEQ_BATCH_SIZE>& b = *(seq[i >> SEQ_BATCH_BIT]);
			i &= SEQ_BATCH_SIZE - 1;
			for (int t = 0; t < 4; t++) {
				cnt[t] -= b[i++];
			}
		}
	};

	template<typename T> struct Array {
		vector<shared_ptr<array<T, SEQ_BATCH_SIZE> > > seq;
		int p = 0;

		Array() : seq(1) {}

		void append(const T& t) {
			array<T, SEQ_BATCH_SIZE> & b = *(seq.back());
			b[p++] = t;
			if (p == SEQ_BATCH_SIZE) {
				p = 0;
				seq.emplace_back();
			}
		}

		const T& get(size_t id) const {
			bitset<SEQ_BATCH_SIZE>& b = *(seq[id >> SEQ_BATCH_BIT]);
			return b[id & (SEQ_BATCH_SIZE - 1)];
		}
	};
	
	struct Gene{
		struct Initializer {
			array<float, 4> pi;
			vector<vector<int> > species2ind;
			int* indSiteRep2kernal;
			size_t* ind2seq;
			float* weight;
			int nInd, nSpecies, nSite, nKernal, nRep;

			Initializer(int nInd, int nSpecies, int nSite, int nKernal, int nRep): species2ind(nSpecies),
					nInd(nInd), nSpecies(nSpecies), nSite(nSite), nKernal(nKernal), nRep(nRep),
					indSiteRep2kernal((nRep == 0) ? nullptr : new int[nInd * nSite * nRep]),
					ind2seq(new size_t[nInd]), weight(new float[nKernal]){}
			
			void setIndSiteRep2kernal(int iInd, int iSite, int iRep, int iKernal){
				indSiteRep2kernal[(iInd * nSite + iSite) * nRep + iRep] = iKernal;
			}

			int* species2indRange() const{
				int* result = new int[nSpecies];
				int cnt = 0;
				for (int i = 0; i < nSpecies; i++){
					cnt += species2ind[i].size();
					result[i] = cnt;
				}
				return result;
			}

			int* indBin() const{
				int cnt = 0;
				for (const vector<int> &arr: species2ind) cnt += arr.size();
				int* result = new int[cnt];
				int i = 0;
				for (const vector<int> &arr: species2ind){
					for (const int &e: arr){
						result[i] = e;
						i++;
					}
				}
				return result;
			}
		};
		
		struct Kernal {
			array<array<unsigned short, 4>, 3> cnt;
			score_t scoreCache = 0;
			float weight = 1;
			bool valid = true;

			void reset(){
				valid = false;
				for (int i = 0; i < 3; i++)
					for (int j = 0; j < 4; j++)
						cnt[i][j] = 0;
			}

			void update(int y, int x, const Sequence &seq, size_t pSeq){
				valid = false;
				if (y != -1) seq.rmv(pSeq, cnt[y]);
				if (x != -1) seq.add(pSeq, cnt[x]);
			}

			score_t score(const array<float, 4> &pi){
				if (valid) return scoreCache;
				valid = true;
				scoreCache = weight * scorePos(cnt, pi);
				return scoreCache;
			}
		};

		const int nInd, nSpecies, nSite, nKernal, nRep;
		const shared_ptr<const int[]> species2indRange, indBin;
		const shared_ptr<const int[]> indSiteRep2kernal;
		const shared_ptr<const size_t[]> ind2seq;
		const shared_ptr<Kernal[]> kernal;
		const array<float, 4> pi;

		Gene(): nInd(0), nSpecies(0), nSite(0), nKernal(0), nRep(0),
				pi({0,0,0,0}), species2indRange(nullptr), indBin(nullptr), 
				indSiteRep2kernal(nullptr), ind2seq(nullptr), kernal(nullptr){}

		Gene(const Initializer& init): nInd(init.nInd), nSpecies(init.nSpecies), nSite(init.nSite), nKernal(init.nKernal), nRep(init.nRep),
				pi(init.pi), species2indRange(init.species2indRange()), indBin(init.indBin()), 
				indSiteRep2kernal(init.indSiteRep2kernal), ind2seq(init.ind2seq), kernal(new Kernal[init.nKernal]) {
			for (int i = 0; i < nKernal; i++) kernal[i].weight = init.weight[i];
		}

		void updateCnt(int i, int y, int x, const Sequence &seq) {
			if (i >= nSpecies) return;
			int indStart = (i == 0) ? 0 : species2indRange[i - 1];
			int indEnd = species2indRange[i];
			for (int indIt = indStart; indIt < indEnd; indIt++) {
				int iInd = indBin[indIt];
				if (nRep == 0){
					for (int iSite = 0; iSite < nSite; iSite++){
						kernal[iSite].update(y, x, seq, ind2seq[iInd] + iSite);
					}
				}
				else {
					for (int iSite = 0; iSite < nSite; iSite++){
						for (int iRep = 0; iRep < nRep; iRep++){
							int iKernal = indSiteRep2kernal[(iInd * nSite + iSite) * nRep + iRep];
							if (iKernal != -1) kernal[iKernal].update(y, x, seq, ind2seq[iInd] + iSite);
						}
					}
				}
			}
		}

		score_t scoreCnt() {
			score_t score = 0;
			for (int iKernal = 0; iKernal < nKernal; iKernal++){
				score += kernal[iKernal].score(pi);
			}
			return score;
		}

		void clearCntScore(){
			for (int iKernal = 0; iKernal < nKernal; iKernal++){
				kernal[iKernal].reset();
			}
		}
	};

	int nThreads = 1, nSpecies = 0;
	vector<Gene> genes;
	Sequence seq;
};

struct Tripartition{
	vector<vector<char> > color;
	TripartitionInitializer& TI;

	Tripartition(TripartitionInitializer &init): TI(init), color(init.nThreads, vector<char>(init.nSpecies, -1)){
		for (TripartitionInitializer::Gene &g: TI.genes){
			g.clearCntScore();
		}
	}

	void updatePart(int part, int x, int i){
		int start = TI.genes.size() * part / TI.nThreads, end = TI.genes.size() * (1 + part) / TI.nThreads;
		int y = color[part][i];
		color[part][i] = x;
		for (int a = start; a < end; a++){
			TI.genes[a].updateCnt(i, y, x, TI.seq);
		}
	}
	
	score_t scorePart(int part){
		int start = TI.genes.size() * part / TI.nThreads, end = TI.genes.size() * (1 + part) / TI.nThreads;
		score_t result = 0;
		for (int a = start; a < end; a++){
			result += TI.genes[a].scoreCnt();
		}
		return result;
	}
};
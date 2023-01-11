#define OBJECTIVE_VERSION "0"

#include<vector>
#include<array>
#include<thread>
#include<future>
#include<bitset>
#include<memory>
#include "threadpool.hpp"

#define CUSTOMIZED_ANNOTATION
#define LOCAL_BOOTSTRAP

using namespace std;

struct CustomizedAnnotation{
	vector<array<score_t,3> > bs;

	CustomizedAnnotation(){}
	CustomizedAnnotation(int len): bs(len){}

	array<int, 3> bootstrap(){
		array<int, 3> result = {};
		for (const array<score_t,3> &e: bs){
			if (e[0] > e[1] && e[0] > e[2]) result[0]++;
			if (e[1] > e[0] && e[1] > e[2]) result[1]++;
			if (e[2] > e[0] && e[2] > e[1]) result[2]++;
		}
		return result;
	}

	CustomizedAnnotation operator+ (const CustomizedAnnotation& o) const{
		CustomizedAnnotation r(bs.size());
		for (int i = 0; i < bs.size(); i++){
			for (int j = 0; j < 3; j++){
				r.bs[i][j] = bs[i][j] + o.bs[i][j];
			}
		}
		return r;
	}

	CustomizedAnnotation& operator+= (const CustomizedAnnotation& o){
		for (int i = 0; i < bs.size(); i++){
			for (int j = 0; j < 3; j++){
				bs[i][j] += o.bs[i][j];
			}
		}
		return *this;
	}
};

inline score_t scorePos(const array<array<unsigned short, 4>, 3> &cnt, score_t pq){
	const unsigned ax0 = cnt[0][0], bx0 = cnt[0][1], xa0 = cnt[0][2], xb0 = cnt[0][3], xx0 = ax0 + bx0;
	const unsigned ax1 = cnt[1][0], bx1 = cnt[1][1], xa1 = cnt[1][2], xb1 = cnt[1][3], xx1 = ax1 + bx1;
	const unsigned ax2 = cnt[2][0], bx2 = cnt[2][1], xa2 = cnt[2][2], xb2 = cnt[2][3], xx2 = ax2 + bx2;
	
	return (xx0 * (xx0 - 1) * pq - ax0 * bx0) * (2 * xx1 * xx2 * pq - (xa1 * xb2 + xb1 * xa2))
		 + (xx0 * (xx0 - 1) * pq - xa0 * xb0) * (2 * xx1 * xx2 * pq - (ax1 * bx2 + bx1 * ax2))
		 + (xx1 * (xx1 - 1) * pq - ax1 * bx1) * (2 * xx2 * xx0 * pq - (xa2 * xb0 + xb2 * xa0))
		 + (xx1 * (xx1 - 1) * pq - xa1 * xb1) * (2 * xx2 * xx0 * pq - (ax2 * bx0 + bx2 * ax0))
		 + (xx2 * (xx2 - 1) * pq - ax2 * bx2) * (2 * xx0 * xx1 * pq - (xa0 * xb1 + xb0 * xa1))
		 + (xx2 * (xx2 - 1) * pq - xa2 * xb2) * (2 * xx0 * xx1 * pq - (ax0 * bx1 + bx0 * ax1));
}

const int SEQ_BATCH_BIT = 25, SEQ_BATCH_SIZE = 1 << SEQ_BATCH_BIT;

struct TripartitionInitializer{
	struct Sequence{
		vector<shared_ptr<bitset<SEQ_BATCH_SIZE> > > seq;
		size_t p = 0;

		Sequence(){
			seq.emplace_back(new bitset<SEQ_BATCH_SIZE>());
		}

		template<bool A, bool C, bool G, bool T> void append(char c1, char c2) {
			bitset<SEQ_BATCH_SIZE>& b = *(seq.back());
			switch (c1) {
				case 'A': case 'a': b.set(p + A, 1); break;
				case 'C': case 'c': b.set(p + C, 1); break;
				case 'G': case 'g': b.set(p + G, 1); break;
				case 'T': case 't': b.set(p + T, 1); break;
			}
			p += 2;
			switch (c2) {
				case 'A': case 'a': b.set(p + A, 1); break;
				case 'C': case 'c': b.set(p + C, 1); break;
				case 'G': case 'g': b.set(p + G, 1); break;
				case 'T': case 't': b.set(p + T, 1); break;
			}
			p += 2;
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
			score_t pq;
			vector<vector<int> > species2ind;
			int* indSiteRep2kernal;
			size_t* ind2seq;
			score_t weight = 1;
			int nInd, nSpecies, nSite, nKernal, nRep;

			Initializer(int nInd, int nSpecies, int nSite, int nKernal, int nRep): species2ind(nSpecies),
					nInd(nInd), nSpecies(nSpecies), nSite(nSite), nKernal(nKernal), nRep(nRep),
					indSiteRep2kernal((nRep == 0) ? nullptr : new int[nInd * nSite * nRep]), ind2seq(new size_t[nInd]){}
			
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

			void reset(){
				for (int i = 0; i < 3; i++)
					for (int j = 0; j < 4; j++)
						cnt[i][j] = 0;
			}

			void update(int y, int x, const Sequence &seq, size_t pSeq){
				if (y != -1) seq.rmv(pSeq, cnt[y]);
				if (x != -1) seq.add(pSeq, cnt[x]);
			}

			score_t score(const score_t pq){
				return scorePos(cnt, pq);
			}
		};

		const int nInd, nSpecies, nSite, nKernal, nRep;
		const shared_ptr<const int[]> species2indRange, indBin;
		const shared_ptr<const int[]> indSiteRep2kernal;
		const shared_ptr<const size_t[]> ind2seq;
		shared_ptr<Kernal[]> kernal;
		const score_t pq;
		const score_t weight = 1;
		
		score_t scoreCache = 0;
		bool valid = false;

		Gene(): nInd(0), nSpecies(0), nSite(0), nKernal(0), nRep(0),
				pq(0), species2indRange(nullptr), indBin(nullptr), 
				indSiteRep2kernal(nullptr), ind2seq(nullptr), kernal(nullptr){}

		Gene(const Initializer& init): nInd(init.nInd), nSpecies(init.nSpecies), nSite(init.nSite), nKernal(init.nKernal), nRep(init.nRep),
				pq(init.pq), species2indRange(init.species2indRange()), indBin(init.indBin()), weight(init.weight),
				indSiteRep2kernal(init.indSiteRep2kernal), ind2seq(init.ind2seq), kernal(new Kernal[init.nKernal]) {}

		void updateCnt(int i, int y, int x, const Sequence &seq) {
			if (i >= nSpecies) return;
			int indStart = (i == 0) ? 0 : species2indRange[i - 1];
			int indEnd = species2indRange[i];
			if (indStart == indEnd) return;
			valid = false;
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
			if (valid) return scoreCache;
			valid = true;
			score_t score = 0;
			for (int iKernal = 0; iKernal < nKernal; iKernal++){
				score += kernal[iKernal].score(pq);
			}
			scoreCache = weight * score;
			return scoreCache;
		}

		void clearCntScore(){
			valid = false;
			for (int iKernal = 0; iKernal < nKernal; iKernal++){
				kernal[iKernal].reset();
			}
		}
	};

	int nThreads = 1, nSpecies = 0;
	vector<Gene> genes;
	Sequence seq;
	int ufb_size = 1000, ufb_fold = 1000;
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

inline score_t quadPos(const array<unsigned short, 4> &cnt0, const array<unsigned short, 4> &cnt1, 
		const array<unsigned short, 4> &cnt2, const array<unsigned short, 4> &cnt3, score_t pq){
	const unsigned ax0 = cnt0[0], bx0 = cnt0[1], xa0 = cnt0[2], xb0 = cnt0[3], xx0 = ax0 + bx0;
	const unsigned ax1 = cnt1[0], bx1 = cnt1[1], xa1 = cnt1[2], xb1 = cnt1[3], xx1 = ax1 + bx1;
	const unsigned ax2 = cnt2[0], bx2 = cnt2[1], xa2 = cnt2[2], xb2 = cnt2[3], xx2 = ax2 + bx2;
	const unsigned ax3 = cnt3[0], bx3 = cnt3[1], xa3 = cnt3[2], xb3 = cnt3[3], xx3 = ax3 + bx3;
	
	return (2 * xx0 * xx1 * pq - (xa0 * xb1 + xb0 * xa1)) * (2 * xx2 * xx3 * pq - (ax2 * bx3 + bx2 * ax3))
		 + (2 * xx0 * xx1 * pq - (ax0 * bx1 + bx0 * ax1)) * (2 * xx2 * xx3 * pq - (xa2 * xb3 + xb2 * xa3));
}

inline array<score_t, 3> quadPos(const array<array<unsigned short, 4>, 4> &cnt, score_t pq){
	return {quadPos(cnt[0], cnt[1], cnt[2], cnt[3], pq),
			quadPos(cnt[0], cnt[2], cnt[1], cnt[3], pq),
			quadPos(cnt[0], cnt[3], cnt[1], cnt[2], pq)};
}

struct Quadrupartition{
	struct Gene{
		struct Kernal {
			array<array<unsigned short, 4>, 4> cnt;
			
			void reset(){
				for (int i = 0; i < 4; i++)
					for (int j = 0; j < 4; j++)
						cnt[i][j] = 0;
			}

			void update(int y, int x, const TripartitionInitializer::Sequence &seq, size_t pSeq){
				if (y != -1) seq.rmv(pSeq, cnt[y]);
				if (x != -1) seq.add(pSeq, cnt[x]);
			}

			array<score_t, 3> score(score_t pq){
				return quadPos(cnt, pq);
			}
		};

		const int nInd, nSpecies, nSite, nKernal, nRep;
		const shared_ptr<const int[]> species2indRange, indBin;
		const shared_ptr<const int[]> indSiteRep2kernal;
		const shared_ptr<const size_t[]> ind2seq;
		const shared_ptr<Kernal[]> kernal;
		const score_t pq;
		const int ufb_size, ufb_fold;
		shared_ptr<int[]> ufb_offset;
		const score_t weight = 1;
		
		array<score_t, 3> scoreCache;
		bool valid = true;

		Gene(TripartitionInitializer::Gene& init, int ufb_size, int ufb_fold): nInd(init.nInd), nSpecies(init.nSpecies), nSite(init.nSite), nKernal(init.nKernal), nRep(init.nRep),
				pq(init.pq), species2indRange(init.species2indRange), indBin(init.indBin), 
				indSiteRep2kernal(init.indSiteRep2kernal), ind2seq(init.ind2seq), kernal(new Kernal[init.nKernal]),
				ufb_size(ufb_size), ufb_fold(ufb_fold), ufb_offset(new int[ufb_fold]), weight(init.weight) {
			for (int i = 0; i < ufb_fold; i++) ufb_offset[i] = rand() % ufb_size;
			init.kernal.reset();
		}

		void updateCnt(int i, int y, int x, const TripartitionInitializer::Sequence &seq) {
			if (i >= nSpecies) return;
			int indStart = (i == 0) ? 0 : species2indRange[i - 1];
			int indEnd = species2indRange[i];
			if (indStart == indEnd) return;
			valid = false;
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

		array<score_t, 3> scoreCnt() {
			valid = true;
			array<score_t, 3> score = {};
			for (int iKernal = 0; iKernal < nKernal; iKernal++){
				array<score_t, 3> temp = kernal[iKernal].score(pq);
				for (int i = 0; i < 3; i++) score[i] += temp[i];
			}
			scoreCache = {score[0] * weight, score[1] * weight, score[2] * weight};
			return scoreCache;
		}

		void annotate(CustomizedAnnotation &annot){
			valid = true;
			array<score_t, 3> score = {};
			vector<array<score_t, 3> > temp(nKernal);
			for (int iKernal = 0; iKernal < nKernal; iKernal++){
				temp[iKernal] = kernal[iKernal].score(pq);
				for (int i = 0; i < 3; i++) score[i] += temp[iKernal][i];
			}
			scoreCache = {score[0] * weight, score[1] * weight, score[2] * weight};
			for (int t = 0; t < ufb_fold; t++){
				int pos = ufb_offset[t];
				for (int iKernal = 0; iKernal < nKernal; iKernal++){
					for (int i = 0; i < 3; i++) annot.bs[pos][i] += temp[iKernal][i];
					pos++;
					if (pos == ufb_size) pos = 0;
				}
			}
		}

		void clearCntScore(){
			valid = false;
			for (int iKernal = 0; iKernal < nKernal; iKernal++){
				kernal[iKernal].reset();
			}
		}
	};

	vector<Gene> genes;
	vector<vector<char> > color;
	TripartitionInitializer& TI;

	Quadrupartition(TripartitionInitializer &init): TI(init), color(init.nThreads, vector<char>(init.nSpecies, -1)){
		for (TripartitionInitializer::Gene &g: TI.genes){
			genes.emplace_back(g, TI.ufb_size, TI.ufb_fold);
			genes.back().clearCntScore();
		}
	}

	void updatePart(int part, int x, int i){
		int start = TI.genes.size() * part / TI.nThreads, end = TI.genes.size() * (1 + part) / TI.nThreads;
		int y = color[part][i];
		color[part][i] = x;
		for (int a = start; a < end; a++){
			genes[a].updateCnt(i, y, x, TI.seq);
		}
	}
	
	void scorePart(int part, array<score_t, 3> &result){
		array<score_t, 3> local_result = {};
		int start = TI.genes.size() * part / TI.nThreads, end = TI.genes.size() * (1 + part) / TI.nThreads;
		for (int a = start; a < end; a++){
			array<score_t, 3> temp = genes[a].scoreCnt();
			for (int i = 0; i < 3; i++) local_result[i] += temp[i];
		}
		result = local_result;
	}
	
	void annotatePart(int part, CustomizedAnnotation &result){
		CustomizedAnnotation local_result(result.bs.size());
		int start = TI.genes.size() * part / TI.nThreads, end = TI.genes.size() * (1 + part) / TI.nThreads;
		for (int a = start; a < end; a++){
			genes[a].annotate(local_result);
		}
		result = local_result;
	}

	void update(int x, int i){
		vector<thread> thrds;
		for (int p = 1; p < color.size(); p++) thrds.emplace_back(&Quadrupartition::updatePart, this, p, x, i);
		updatePart(0, x, i);
		for (thread &t: thrds) t.join();
	}
	
	array<double, 3> score(){
		array<double, 3> res;
		vector<array<score_t, 3> > t(color.size());
		vector<thread> thrds;
		for (int p = 1; p < color.size(); p++) thrds.emplace_back(&Quadrupartition::scorePart, this, p, ref(t[p]));
		scorePart(0, t[0]);
		for (thread &t: thrds) t.join();
		res[0] = 0;
		res[1] = 0;
		res[2] = 0;
		for (int p = 0; p < color.size(); p++){
			res[0] += t[p][0];
			res[1] += t[p][1];
			res[2] += t[p][2];
		}
		return res;
	}

	CustomizedAnnotation annotate(){
		CustomizedAnnotation res(TI.ufb_size);
		vector<CustomizedAnnotation> t(color.size(), TI.ufb_size);
		vector<thread> thrds;
		for (int p = 1; p < color.size(); p++) thrds.emplace_back(&Quadrupartition::annotatePart, this, p, ref(t[p]));
		annotatePart(0, t[0]);
		for (thread &t: thrds) t.join();
		for (int p = 0; p < color.size(); p++){
			res += t[p];
		}
		return res;
	}
};
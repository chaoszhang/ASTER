#define DRIVER_VERSION "5"

/* CHANGE LOG
 * 5: Early termination for high coverage read files
 * 4: Parallelization
 * 3: A different strategy
 * 2: Final score normalization
 * 1: Remove non-effective sites
 */

#include<iostream>
#include<fstream>
#include<unordered_map>
#include<unordered_set>
#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<cstring>
#include<algorithm>
#include<random>
#include<thread>
#include<mutex>

//#define DEBUGINFO

#define ROOTING
#define NAME_MAPPING

//#define CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH
#ifdef CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH
//#define CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH_ROUGH
#else
#define CUSTOMIZED_ANNOTATION_LENGTH
#endif

//#define LARGE_DATA
#ifdef LARGE_DATA
typedef long double score_t;
typedef long long count_t;
#else
typedef double score_t;
typedef int count_t;
#endif

#include "sequence.hpp"
#include "algorithms.hpp"
#include "sequtils.hpp"

using namespace std;

const char BITS2CHAR[4] = {'A', 'C', 'G', 'T'};

unsigned char CHAR2BITS(char c){
    if (c == 'A') return 0;
    if (c == 'C') return 1;
    if (c == 'G') return 2;
    if (c == 'T') return 3;
    return 0;
}

unsigned char CHAR2MASK(char c){
    if (c == 'N') return 1;
    return 0;
}

char REVERSE_COMPLEMENT(char c){
    if (c == 'A') return 'T';
    if (c == 'C') return 'G';
    if (c == 'G') return 'C';
    if (c == 'T') return 'A';
    return 'N';
}

template<size_t K> size_t REVERSE_COMPLEMENT(size_t b){
    size_t res = 0;
    for (size_t i = 0; i < K; i++){
        res <<= 2;
        res |= (3 ^ (b & 3));
        b >>= 2;
    }
    return res;
}

template<size_t K> string print(size_t b){
    string res;
    for (int i = K - 1; i >= 0; i--){
        res += BITS2CHAR[(b >> (2 * i)) & 3];
    }
    return res;
}

template<size_t K> string print(size_t b1, size_t b2){
    return print<K>(b1) + "-" + SeqParser::STRING_REVERSE_COMPLEMENT(print<K>(b2));
}

template<size_t K> string print(size_t b1, char c, size_t b2){
    return print<K>(b1) + c + SeqParser::STRING_REVERSE_COMPLEMENT(print<K>(b2));
}

template<size_t K> class KmerTable{
    constexpr static const size_t MASK = (1LL << (K * 2)) - 1;

    uint8_t* table;
    vector<vector<size_t> > cache;
    size_t cacheSize = 0, nThreads;

public:
    constexpr static const size_t LEN = (1LL << (K * 4 - 1)) - (1LL << (K * 2 - 1));
    constexpr static const size_t CACHE_MAX = (1LL << 27);

    size_t fillcnt = 0;

    static void shift(size_t &bits, size_t &mask, const unsigned char c){
        bits = ((bits << 2) | CHAR2BITS(c)) & MASK;
        mask = ((mask << 2) | CHAR2MASK(c)) & MASK;
    }

    static void rshift(size_t &bits, size_t &mask, const unsigned char c){
        bits = (bits >> 2) | ((3LL ^ CHAR2BITS(c)) << (K * 2 - 2));
        mask = ((mask << 2) | CHAR2MASK(c)) & MASK;
    }

    static size_t REVERSE_COMPLEMENT(size_t bits){
        size_t res = 0;
        for (size_t i = 0; i < 2 * K; i++){
            res = (res << 2) | (3 ^ (bits & 3));
            bits >>= 2;
        }
        return res;
    }

    KmerTable(size_t nThreads = 1): table(new uint8_t[LEN]{}), cache(nThreads), nThreads(nThreads){}

    ~KmerTable(){
        delete table;
    }

    /*
    void addChar(size_t p, uint8_t c) {
        uint8_t b = (table[p] & 7);
        if (b == 3 || b == 4 + c) return;
        if (b == 0) table[p] += 4 + c;
        else table[p] += 3 - b;
    }
    */

    void addChar(size_t p, uint8_t c) {
        size_t threadID = (p * nThreads) / LEN;
        cache[threadID].push_back(p * 4 + c);
        cacheSize++;
        if (cacheSize >= CACHE_MAX) performAddChar();
    }

    void performAddChar(){
        vector<thread> threads;
        for (size_t i = 1; i < nThreads; i++){
            threads.emplace_back(&KmerTable<K>::performAddCharThread, this, i);
        }
        performAddCharThread(0);
        for (thread &t: threads) t.join();
        for (size_t i = 0; i < nThreads; i++){
            cache[i].clear();
        }
        cacheSize = 0;
    }


    /* bit table (last 3 bits) for central-mer:
     * 0 - vacancy
     * 3 - mutiple hits (error)
     * 4 - A
     * 5 - C
     * 6 - G
     * 7 - T
     */
    void performAddCharThread(size_t threadID){
        for (size_t e: cache[threadID]){
            size_t p = e / 4;
            uint8_t c = (e & 3);
            uint8_t b = (table[p] & 7);
            if (b == 3 || b == 4 + c) continue;
            if (b == 0) table[p] += 4 + c;
            else table[p] += 3 - b;
        }
    }

    void add(const string &seqs, const string &seqn, bool bothStrands = true){
        size_t b1 = 0, b2 = 0;
        size_t m1 = 0, m2 = 0;
        for (size_t i = 0; i + K + 1 < seqn.size(); i++){
            shift(b1, m1, seqn[i]);
            rshift(b2, m2, seqn[i + K + 1]);
            char c = seqs[i + 1];
            if (i >= K - 1 && CHAR2MASK(c) == 0){
                if (b1 > b2 && m1 == 0 && m2 == 0){
                    size_t p = b1 * (b1 - 1) / 2 + b2;
                    addChar(p, CHAR2BITS(c));
                }
            }
        }
        if (bothStrands) add(SeqParser::STRING_REVERSE_COMPLEMENT(seqs), SeqParser::STRING_REVERSE_COMPLEMENT(seqn), false);
    }

    void postprocess(){
        performAddChar();
        fillcnt = 0;
        vector<size_t> cnts(nThreads);
        vector<thread> threads;
        for (size_t i = 1; i < nThreads; i++) {
            threads.emplace_back(&KmerTable<K>::performPostprocess, this, i, ref(cnts[i]));
        }
        performPostprocess(0, cnts[0]);
        for (thread& t : threads) t.join();
        for (size_t i = 0; i < nThreads; i++) {
            fillcnt += cnts[i];
        }
    }

    /* bit table (multiple of 8):
     * 8 - mutiple hits (error) in previous taxa
     * multiple of 16: valid counts in previous taxa
     */
    void performPostprocess(size_t threadID, size_t &realcnt) {
        size_t cnt = 0;
        size_t start = (threadID * LEN) / nThreads;
        size_t end = ((threadID + 1) * LEN) / nThreads;
        for (size_t i = start; i < end; i++) {
            if (table[i] == 0) continue;
            cnt++;
            if ((table[i] & 7) == 0) continue;
            if ((table[i] & 7) == 3) {
                if ((table[i] & 8) == 0) table[i] += 8 - 3;
            }
            else {
                table[i] &= 248;
                if (table[i] < 240) table[i] += 16;
            }
        }
        realcnt = cnt;
    }

    double fillProportion() const{
        return fillcnt / (double) LEN;
    }

    vector<size_t> frequentPatterns(size_t n) const{
        size_t freq[16] = {};
        for (size_t i = 0; i < LEN; i++){
            if (table[i] == 0) continue;
            if ((table[i] & 15) == 11) continue;
            size_t level = table[i] / 16;
            freq[level]++;
        }

        size_t threshold, num = n, denom = freq[0];
        for (threshold = 15; threshold > 0; threshold--) {
            if (num <= freq[threshold]) {
                denom = freq[threshold];
                break;
            }
            num -= freq[threshold];
        }

        vector<size_t> res;
        for (size_t i = 0; i < LEN; i++) {
            if (table[i] == 0) continue;
            if ((table[i] & 15) == 11) continue;
            size_t level = table[i] / 16;
            if (level > threshold) res.push_back(i);
            if (level == threshold) {
                if (denom == 0) break;
                if (i % denom < num) {
                    res.push_back(i);
                    num--;
                }
                denom--;
            }
        }

        LOG << "Threshold: " << threshold << endl;
        
        #ifdef DEBUGINFO
        for (size_t e: res){
            size_t a1 = (static_cast<size_t>(floor(sqrt(8 * e + 1))) + 1) / 2;
            size_t a2 = e - a1 * (a1 - 1) / 2;
            cerr << print<K>(a1, a2) << " ";
        }
        cerr << endl;
        #endif

        return res;
    }
};

template<size_t K> class KmerConsensus {
    constexpr static const size_t SHIFT = 8;
    constexpr static const size_t M = (1LL << (K * 4 - 1 - SHIFT)) - (1LL << (K * 2 - 1 - SHIFT));
    constexpr static const size_t MASK = (1LL << (K * 2)) - 1;

    constexpr static const size_t LEN = (1LL << (K * 4 - 1)) - (1LL << (K * 2 - 1));
    constexpr static const size_t CACHE_MAX = (1LL << 26);

    struct Profile {
        int32_t sample[7] = { -1, -1, -1, -1, -1, -1, -1 };
        int16_t cnt = 0;
        bool sampled = false;

        void add(size_t b1, size_t b2) {
            if (sampled) return;
            sampled = true;
            cnt++;
            b1 >>= 2 * K;
            b2 >>= 2 * K;
            int32_t b = (b1 << (2 * (16 - K))) + b2;
            for (int i = 0; i < 7; i++) {
                if (sample[i] == -1) {
                    sample[i] = b;
                    break;
                }
            }
        }

        pair<size_t, double> output(size_t p) {
            size_t a1 = (static_cast<int>(floor(sqrt(8 * p + 1))) + 1) / 2;
            size_t a2 = p - a1 * (a1 - 1) / 2;
            size_t b = 0;
            if (a2 >= a1) {
                cerr << "What? a1 <= a2!\n";
                exit(-1);
            }
            int nSample = 0, totalCnt = 0;
            for (int j = 0; j < 7; j++) {
                if (sample[j] == -1) break;
                nSample++;
            }
            if (nSample == 0) {
                return { -1, 0.0 };
            }
            for (int i = 31 - K * 2; i >= 0; i--) {
                int cntC[4] = {}, bestC = 0;
                for (int j = 0; j < nSample; j++) {
                    cntC[(sample[j] >> (2 * i)) & 3]++;
                }
                for (int j = 1; j < 4; j++) {
                    if (cntC[j] > cntC[bestC]) bestC = j;
                }
                b = 4 * b + bestC;
                totalCnt += cntC[bestC];
            }
            double temp = totalCnt / (16.0 - K) / nSample - 1;
            if (temp < 0) return { -1, 0.0 };
            size_t b1 = (b >> (32 - 2 * K)), b2 = b - (b1 << (32 - 2 * K));
            size_t ba1 = (b1 << (2 * K)) + a1;
            size_t ba2 = (b2 << (2 * K)) + a2;
            if (ba1 > ba2) return { (ba1 << 32) + ba2, temp * cnt};
            else return { (ba2 << 32) + ba1, temp * cnt };
        }
    };

public:
    vector<size_t> sortedPattern, pattern2pos;
    vector<Profile> profile;
    vector<vector<array<size_t, 3> > > cache;
    size_t cacheSize = 0, nThreads;

    KmerConsensus(const vector<size_t> &sortedPattern, size_t nThreads = 1): 
            sortedPattern(sortedPattern), profile(sortedPattern.size()), cache(nThreads), nThreads(nThreads) {
        pattern2pos.reserve(M + 1);
        size_t n = sortedPattern.size();
        for (size_t i = 0, j = 0; i < M + 1; i++) {
            while (j < n && sortedPattern[j] < (i << SHIFT)) j++;
            pattern2pos.push_back(j);
        }
    }
    
    void profileAdd(size_t pattern, size_t b1, size_t b2){
        size_t threadID = (pattern * nThreads) / LEN;
        cache[threadID].push_back({pattern, b1, b2});
        cacheSize++;
        if (cacheSize >= CACHE_MAX) performProfileAdd();
    }

    void performProfileAdd(){
        vector<thread> threads;
        for (size_t i = 1; i < nThreads; i++){
            threads.emplace_back(&KmerConsensus<K>::performProfileAddThread, this, i);
        }
        performProfileAddThread(0);
        for (thread &t: threads) t.join();
        for (size_t i = 0; i < nThreads; i++){
            cache[i].clear();
        }
        cacheSize = 0;
    }

    void performProfileAddThread(size_t threadID){
        for (const array<size_t, 3> &e: cache[threadID]){
            size_t pattern = e[0];
            size_t L = pattern2pos[pattern >> SHIFT], R = pattern2pos[(pattern >> SHIFT) + 1];
            size_t i = L;
            while (i < R && sortedPattern[i] < pattern) i++;
            if (sortedPattern[i] == pattern) profile[i].add(e[1], e[2]);
        }
    }

    /*
    size_t locate(size_t pattern) {
        size_t L = pattern2pos[pattern >> SHIFT], R = pattern2pos[(pattern >> SHIFT) + 1];
        size_t i = L;
        while (i < R && sortedPattern[i] < pattern) i++;
        if (sortedPattern[i] == pattern) return i;
        return -1;
    }
    */
    
    void add(const string& seqs, const string& seqn, bool bothStrands = true) {
        size_t b1 = 0, b2 = 0;
        size_t m1 = 0, m2 = 0;
        for (size_t i = 0; i + 17 < seqs.size(); i++) {
            KmerTable<16>::shift(b1, m1, seqn[i]);
            KmerTable<16>::rshift(b2, m2, seqn[i + 17]);
            char c = seqs[i + 1];
            if (i >= 15 && CHAR2MASK(c) == 0) {
                size_t mb1 = (b1 & MASK);
                size_t mb2 = (b2 & MASK);
                if (mb1 > mb2 && m1 == 0 && m2 == 0) {
                    size_t p = mb1 * (mb1 - 1) / 2 + mb2;
                    profileAdd(p, b1, b2);
                }
            }
        }
        if (bothStrands) add(SeqParser::STRING_REVERSE_COMPLEMENT(seqs), SeqParser::STRING_REVERSE_COMPLEMENT(seqn), false);
    }

    void postprocess() {
        performProfileAdd();
        for (size_t i = 0; i < sortedPattern.size(); i++) {
            if (profile[i].sampled) profile[i].sampled = false;
        }
    }

    vector<pair<size_t, double> > generateConsensus(size_t n) {
        vector<pair<size_t, double> > temp, res;
        for (size_t i = 0; i < sortedPattern.size(); i++) {
            pair<size_t, double> e = profile[i].output(sortedPattern[i]);
            if (e.second > 0) temp.push_back(e);
        }
        sort(temp.begin(), temp.end(), [](const pair<size_t, double> &a, const pair<size_t, double> &b){return a.second > b.second;});
        for (size_t i = 0; i < n && i < temp.size(); i++){
            res.push_back(temp[i]);
        }

        #ifdef DEBUGINFO
        for (const auto& e: res){
            cerr << print<16>(e.first >> 32, e.first & 0xffffffff) << ":" << e.second << " ";
        }
        cerr << endl;
        #endif

        return res;
    }
};

class SNP {
public:
    struct SEQ{
        string seq;
        size_t uniqCnt = 0;

        #ifdef DEBUGINFO
        vector<size_t> debug;

        void debuginfo(){
            for (int i = 0; i < seq.size(); i++){
                if (seq[i] == '\0' || seq[i] == 'N') cerr << "----------------N---------------- ";
                else cerr << print<16>(debug[i] >> 32, seq[i], debug[i] & 0xffffffff) << " ";
            }
            cerr << endl;
        }
        #endif

        SEQ(){}

        SEQ(int n){
            for (int i = 0; i < n; i++){
                seq += '\0';

                #ifdef DEBUGINFO
                debug.push_back(-1);
                #endif
            }
        }

        void set(int i, char c){
            if (seq[i] == c || seq[i] == 'N') return;
            uniqCnt++;
            if (seq[i] == '\0') seq[i] = c;
            else seq[i] = 'N';
        }

        string get(){
            for (size_t i = 0; i < seq.length(); i++){
                if (seq[i] == '\0') seq[i] = 'N';
            }
            return seq;
        }

        size_t uniqueCount(){
            return uniqCnt;
        }

        size_t countAT(){
            size_t cnt = 0;
            for (size_t i = 0; i < seq.length(); i++){
                if (seq[i] == 'A' || seq[i] == 'T') cnt++;
            }
            return cnt;
        }

        size_t countCG(){
            size_t cnt = 0;
            for (size_t i = 0; i < seq.length(); i++){
                if (seq[i] == 'C' || seq[i] == 'G') cnt++;
            }
            return cnt;
        }
    };

private:
    constexpr static const size_t MASK = 0xffffffff;
    constexpr static const size_t NMASKS = 14, MAXDIST = 3;
    const size_t masks[NMASKS] = {0x33333333, 0xcccccccc, 0x0f0f0f0f, 0xf0f0f0f0, 0x00ff00ff, 0xff00ff00, 0x0000ffff, 0xffff0000,
                                  0x3333cccc, 0xcccc3333, 0x0f0ff0f0, 0xf0f00f0f, 0x00ffff00, 0xff0000ff};
    vector<size_t> patterns;
    vector<double> weights;
    vector<int> unionfind;
    unordered_map<size_t, int> masked[NMASKS];
    vector<int> id2pos;
    int npos = 0;
    
    static size_t dist(size_t a, size_t b) {
        size_t res = 0;
        size_t d = (a ^ b);
        for (size_t i = 0; i < 32; i++) {
            if (d & 3) res++;
            d >>= 2;
        }
        return res;
    }

    int ufroot(int node) {
        if (unionfind[node] == -1) return node;
        return unionfind[node] = ufroot(unionfind[node]);
    }

    int ufrootConst(int node) const{
        if (unionfind[node] == -1) return node;
        return ufrootConst(unionfind[node]);
    }

    static size_t pairup(size_t b1, size_t b2) {
        if (b1 > b2) return (b1 << 32) + b2;
        else return (b2 << 32) + b1;
    }

public:

    SNP(const vector<pair<size_t, double> > &consensus) {
        for (auto& e : consensus) {
            patterns.push_back(e.first);
            weights.push_back(e.second);
            unionfind.push_back(-1);
        }
        for (size_t i = 0; i < patterns.size(); i++) {
            size_t b1 = (patterns[i] >> 32);
            size_t b2 = (patterns[i] & MASK);
            for (size_t k = 0; k < NMASKS; k++) {
                size_t mb1 = (masks[k] & b1);
                size_t mb2 = (masks[k] & b2);
                size_t mb = pairup(mb1, mb2);
                if (masked[k].count(mb)) {
                    size_t j = masked[k][mb];
                    if (dist(patterns[i], patterns[j]) <= MAXDIST || dist(patterns[i], REVERSE_COMPLEMENT<32>(patterns[j])) <= MAXDIST 
                            || dist(patterns[ufroot(i)], patterns[ufroot(j)]) <= MAXDIST || dist(patterns[ufroot(i)], REVERSE_COMPLEMENT<32>(patterns[ufroot(j)])) <= MAXDIST) {
                        if (weights[ufroot(i)] > weights[ufroot(j)]) {
                            unionfind[ufroot(j)] = ufroot(i);
                            masked[k][mb] = i;
                        }
                        else if (ufroot(i) != ufroot(j)) {
                            unionfind[ufroot(i)] = ufroot(j);
                        }
                    }
                    else {
                        if (weights[ufroot(i)] > weights[ufroot(j)]) {
                            masked[k][mb] = i;
                        }
                    }
                }
                else {
                    masked[k][mb] = i;
                }
            }
        }
        for (size_t i = 0; i < patterns.size(); i++) {
            if (unionfind[i] == -1) id2pos.push_back(npos++);
            else id2pos.push_back(-1);
        }

        #ifdef DEBUGINFO
        for (size_t i = 0; i < patterns.size(); i++) {
            if (unionfind[i] == -1) cerr << print<16>(patterns[i] >> 32, patterns[i] & 0xffffffff) << " ";
        }
        cerr << endl;
        #endif
    }

    SNP(const string &file){
        ifstream fin(file);
        size_t LEN;
        for (size_t k = 0; k < NMASKS; k++) {
            fin >> LEN;
            for (size_t i = 0; i < LEN; i++){
                size_t first;
                size_t second;
                fin >> first >> second;
                masked[k][first] = second;
            }
        }
        fin >> LEN;
        for (size_t i = 0; i < LEN; i++){
            size_t e;
            fin >> e;
            patterns.push_back(e);
        }
        fin >> LEN;
        for (size_t i = 0; i < LEN; i++){
            double e;
            fin >> e;
            weights.push_back(e);
        }
        fin >> LEN;
        for (size_t i = 0; i < LEN; i++){
            int e;
            fin >> e;
            unionfind.push_back(e);
        }
        fin >> LEN;
        for (size_t i = 0; i < LEN; i++){
            size_t e;
            fin >> e;
            id2pos.push_back(e);
        }
        fin >> npos;
    }
    
    void save(const string &file){
        ofstream fout(file);
        for (size_t k = 0; k < NMASKS; k++) {
            fout << masked[k].size() << endl;
            for (const auto &e: masked[k]) fout << e.first << " " << e.second << endl;
        }
        fout << patterns.size() << endl;
        for (size_t e: patterns) fout << e << endl;
        fout << weights.size() << endl;
        for (double e: weights) fout << e << endl;
        fout << unionfind.size() << endl;
        for (int e: unionfind) fout << e << endl;
        fout << id2pos.size() << endl;
        for (int e: id2pos) fout << e << endl;
        fout << npos << endl;
    }

    SEQ initSeq(){
        SEQ seq(npos);
        return seq;
    }

    void add(SEQ& seq, const string &seqs, const string &seqn) const{
        size_t b1 = 0, b2 = 0;
        size_t m1 = 0, m2 = 0;
        for (size_t i = 0; i + 17 < seqs.size(); i++){
            KmerTable<16>::shift(b1, m1, seqn[i]);
            KmerTable<16>::rshift(b2, m2, seqn[i + 17]);
            if (i < 15 || CHAR2MASK(seqs[i + 1]) != 0) continue;
            if (m1 != 0 || m2 != 0) continue;
            size_t b = (b1 << 32) + b2, rb = (b2 << 32) + b1;
            int best = -1;
            char c;
            for (size_t k = 0; k < NMASKS; k++) {
                size_t mb1 = (masks[k] & b1);
                size_t mb2 = (masks[k] & b2);
                size_t mb = pairup(mb1, mb2);
                if (masked[k].count(mb) == 0) continue;
                int temp = masked[k].at(mb), next = ufrootConst(temp);
                if ((dist(patterns[temp], b) <= MAXDIST || dist(patterns[next], b) <= MAXDIST) && (best != -1 || weights[best] < weights[next])) {
                    best = next;
                    c = seqs[i + 1];
                }
                if ((dist(patterns[temp], rb) <= MAXDIST || dist(patterns[next], rb) <= MAXDIST) && (best != -1 || weights[best] < weights[next])) {
                    best = next;
                    c = REVERSE_COMPLEMENT(seqs[i + 1]);
                }
            }
            if (best == -1) continue;
            seq.set(id2pos[best], c);

            #ifdef DEBUGINFO
            seq.debug[id2pos[best]] = (c != seqs[i + 1]) ? rb : b;
            #endif
        }
    }
};

struct Workflow {
    MetaAlgorithm meta;
    TripartitionInitializer& tripInit = meta.tripInit;

    vector<string>& names = meta.names;
    unordered_map<string, int>& name2id = meta.name2id;
    
    int intqcs, intqcn;
    bool eofFlag = false;

    double GCcontent = 0.5;
    size_t nSNP = 0;
    vector<int> ind2species;

    void addName(const string& name) {
        if (name2id.count(name) == 0) {
            name2id[name] = names.size();
            names.push_back(name);
        }
    }

    void buildGeneSeq(TripartitionInitializer::Gene::Initializer &gene, size_t pos, size_t nSite, size_t offset){
        int nInd = ind2species.size();
        array<size_t, 4> cnt = {};
        for (int iInd = 0; iInd < nInd; iInd++) {
            size_t pSeq = pos + iInd * offset;
            for (size_t iSite = 0; iSite < nSite; iSite++) {
                tripInit.seq.add(pSeq + iSite, cnt);
            }
            gene.species2ind[ind2species[iInd]].push_back(iInd);
            gene.ind2seq[iInd] = pSeq;
        }
        gene.pi[0] = (1 - GCcontent) * 0.5;
        gene.pi[1] = GCcontent * 0.5;
        gene.pi[2] = GCcontent * 0.5;
        gene.pi[3] = (1 - GCcontent) * 0.5;
    }

    void formatGene(size_t pos, size_t nSite, size_t offset) {
        size_t nInd = ind2species.size(), nSpecies = names.size(), nKernal = nSite, nRep = 0;
        TripartitionInitializer::Gene::Initializer gene(nInd, nSpecies, nSite, nKernal, nRep);
        buildGeneSeq(gene, pos, nSite, offset);
        tripInit.genes.emplace_back(gene);
    }

    static int countAT(const string seq){
        int cnt = 0;
        for (const char c: seq){
            if (c == 'A' || c == 'T') cnt++;
        }
        return cnt;
    }

    static int countCG(const string seq){
        int cnt = 0;
        for (const char c: seq){
            if (c == 'C' || c == 'G') cnt++;
        }
        return cnt;
    }
    
    void readAlignment(){
        size_t freqAT = 0, freqCG = 0, npos = 0;
        string line;
        vector<int8_t> cntAG, cntCT;
        {
            ifstream fin(ARG.getStringArg("input"));
            while (getline(fin, line)){
                if (line[0] == '>'){
                    string name = line.substr(1);
                    LOG << "Processing " << name << " ...\n";
                    continue;
                }
                if (cntAG.size() == 0){
                    cntAG.resize(line.size());
                    cntCT.resize(line.size());
                }
                freqAT += countAT(line);
                freqCG += countCG(line);
                for (size_t i = 0; i < line.size(); i++){
                    if ((line[i] == 'A' || line[i] == 'G') && cntAG[i] < 2) cntAG[i]++;
                    if ((line[i] == 'C' || line[i] == 'T') && cntCT[i] < 2) cntCT[i]++; 
                }
            }
        }
        {
            ifstream fin(ARG.getStringArg("input"));
            while (getline(fin, line)){
                if (line[0] == '>'){
                    string name = line.substr(1);
                    LOG << "SNP calling " << name << " ...\n";
                    addName(meta.mappedname(name));
                    ind2species.push_back(name2id[meta.mappedname(name)]);
                    continue;
                }
                string filtered;
                for (size_t i = 0; i < line.size(); i++){
                    #ifdef CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH
                    filtered.push_back(line[i]);
                    #else
                    if (cntAG[i] >= 2 && cntCT[i] >= 2) filtered.push_back(line[i]);
                    #endif
                }
                npos += filtered.size();
                for (const char c: filtered) tripInit.seq.append(c);
            }
        }
        GCcontent = ((double) freqCG) / (freqAT + freqCG);
        nSNP = npos / ind2species.size();
    }

    Workflow(int argc, char** argv){
        //string mappingFile;
        meta.initialize(argc, argv);
        LOG << "Make sure you have run 'waster-site -h', read about '-k' command, and ensured you have enough memory to proceed!\n";
        if (ARG.getStringArg("root") != "") addName(ARG.getStringArg("root"));
        tripInit.nThreads = meta.nThreads;

        if (ARG.getIntArg("mode") <= 3){
            intqcs = ARG.getIntArg("qcs");
            if (intqcs >= 94){
                cerr << "Bad quality control threshold for SNPs!\n";
                exit(0);
            }
            LOG << "Quality control: Masking all SNP bases with quality lower than '" << SeqUtils::QUALITY2ASCII[intqcs] << "' for FASTQ inputs.\n";
            intqcn = ARG.getIntArg("qcn");
            if (intqcn >= 94){
                cerr << "Bad quality control threshold for non-SNPs!\n";
                exit(0);
            }
            LOG << "Quality control: Masking all non-SNP bases with quality lower than '" << SeqUtils::QUALITY2ASCII[intqcn] << "' for FASTQ inputs.\n";
            
            #ifdef DEBUGINFO
            init<8>();
            #else
            init<9>();
            #endif
        }
        readAlignment();

        size_t nChunk = meta.nThreads;
        size_t pos = 0, len = nSNP;
        for (size_t i = 0; i < nChunk; i++) {
            size_t s = i * len / nChunk, t = (i + 1) * len / nChunk;
            formatGene(pos + s, t - s, len);
        }

        tripInit.nSpecies = names.size();
    }

    static void processingFile(const SNP& snp, SNP::SEQ& snpseq, string file, int intqcs, int intqcn){
        string log = string("Processing ") + file + " ...\n";
        size_t lastCnt = 0, curCnt = 0, totalLen = 0, nextThreshold = (1LL << 30);
        double thresholdRatio = 1.01;
        LOG << log;
        SeqParser seq(file);
        while (seq.nextSeq()){
            string seqs = seq.getSeq(intqcs), seqn = seq.getSeq(intqcn);
            snp.add(snpseq, seqs, seqn);
            totalLen += seqs.size();
            if (totalLen > nextThreshold){
                nextThreshold += (1LL << 30);
                curCnt = snpseq.uniqueCount();
                if (lastCnt * thresholdRatio > curCnt){
                    LOG << string("Early termination of ") + file + " after " + to_string(totalLen) + " base-pairs as the coverage level is sufficiently high!\n";
                    return;
                }
                lastCnt = curCnt;
            }
        }
    }
    
    template<size_t K> void init(){
        vector<string> indNames, files;
        {
            ifstream fin(ARG.getStringArg("input"));
            string name, file;
            while (fin >> name){
                indNames.push_back(name);
                addName(meta.mappedname(name));
                fin >> file;
                files.push_back(file);
            }
            if (files.size() == 0){
                cerr << "Error: Input list file '" << ARG.getStringArg("input") <<"' empty! Typo?\n";
                exit(0);
            }
        }
        SNP *pSNP;
        if (ARG.getStringArg("continue") == ""){
            size_t nSample = ARG.getIntArg("sampled");
            vector<pair<size_t, double> > generatedConsensus;
            unordered_set<size_t> selected;
            vector<size_t> frequentPatterns;

            {
                KmerTable<K> table(ARG.getIntArg("thread"));
                for (size_t i = 0; i < files.size() && i < nSample; i++){
                    size_t cur = rand() % files.size();
                    while (selected.count(cur)) cur = rand() % files.size();
                    selected.insert(cur);
                    LOG << "Species " << indNames[cur] << " is selected to count the most frequent patterns.\n";
                    SeqParser seq(files[cur]);
                    while (seq.nextSeq()){
                        string seqs = seq.getSeq(intqcs), seqn = seq.getSeq(intqcn);
                        table.add(seqs, seqn);
                    }
                    table.postprocess();
                }
                #ifdef DEBUGINFO
                frequentPatterns = table.frequentPatterns(1000);
                #else
                frequentPatterns = table.frequentPatterns(ARG.getIntArg("pattern"));
                #endif
            }

            {
                selected.clear();
                KmerConsensus<K> consensus(frequentPatterns, ARG.getIntArg("thread"));
                for (size_t i = 0; i < files.size() && i < nSample; i++){
                    size_t cur = rand() % files.size();
                    while (selected.count(cur)) cur = rand() % files.size();
                    selected.insert(cur);
                    LOG << "Species " << indNames[cur] << " is selected to count to build consensus.\n";
                    SeqParser seq(files[cur]);
                    while (seq.nextSeq()){
                        string seqs = seq.getSeq(intqcs), seqn = seq.getSeq(intqcn);
                        consensus.add(seqs, seqn);
                    }
                    consensus.postprocess();
                }
                frequentPatterns.clear();

                #ifdef DEBUGINFO
                generatedConsensus = consensus.generateConsensus(100);
                #else
                generatedConsensus = consensus.generateConsensus(ARG.getIntArg("consensus"));
                #endif
            }

            LOG << generatedConsensus.size() << " consensus sequences selected.\n";

            pSNP = new SNP(generatedConsensus);

            if (ARG.getIntArg("mode") == 2){
                if (ARG.getStringArg("output") == "<standard output>"){
                    cerr << "You cannot save pattern without -o option!\n";
                }
                else {
                    pSNP->save(ARG.getStringArg("output"));
                }
                exit(0);
            }
        }
        else {
            pSNP = new SNP(ARG.getStringArg("continue"));
        }
        SNP &snp = *pSNP;
        string snpName = (ARG.getStringArg("output") != "<standard output>") ? ARG.getStringArg("output") + ".snps.fa" : "snps.fa";
        if (ARG.getIntArg("mode") == 3){
            if (ARG.getStringArg("output") != "<standard output>") ofstream fout(ARG.getStringArg("output"));
        }
        if (ARG.getIntArg("mode") < 3){
            ofstream fout(snpName);
        }
        
        {
            size_t freqAT = 0, freqCG = 0;
            size_t nThreads = ARG.getIntArg("thread");
            for (size_t batch = 0; batch < files.size(); batch += nThreads){
                vector<thread> threads;
                vector<SNP::SEQ> snpseqs;
                for (size_t i = batch; i < files.size() && i < batch + nThreads; i++){
                    snpseqs.push_back(snp.initSeq());
                }
                for (size_t i = batch + 1; i < files.size() && i < batch + nThreads; i++){
                    threads.emplace_back(Workflow::processingFile, ref(snp), ref(snpseqs[i - batch]), files[i], intqcs, intqcn);
                }
                processingFile(snp, snpseqs[0], files[batch], intqcs, intqcn);
                for (size_t i = 0; i < threads.size(); i++){
                    threads[i].join();
                }
                for (size_t i = 0; i < snpseqs.size(); i++){
                    string s = snpseqs[i].get();
                    nSNP = s.size();
                    freqAT += snpseqs[i].countAT();
                    freqCG += snpseqs[i].countCG();
                    string text = string(">") + indNames[batch + i] + "\n" + s + "\n";
                    if (ARG.getIntArg("mode") == 3){
                        if (ARG.getStringArg("output") == "<standard output>"){
                            cout << text;
                        }
                        else {
                            ofstream fout(ARG.getStringArg("output"), ios_base::app);
                            fout << text;
                        }
                    }
                    else {
                        ofstream fout(snpName, ios_base::app);
                        fout << text;
                    }

                    #ifdef DEBUGINFO
                    snpseqs[i].debuginfo();
                    #endif
                }
            }
            GCcontent = ((double) freqCG) / (freqAT + freqCG);
            cerr << "%GC = " << GCcontent << endl;
            if (ARG.getIntArg("mode") == 3) exit(0);

            ARG.getStringArg("input") = snpName;
        }

        /*
        SNP<K>::staticFilter();
        nSNP = freqPatterns.size();
        LOG << nSNP << " effective SNPs found.\n";
        {
            SNP<K> snp;
            for (size_t i = 0; i < files.size(); i++){
                LOG << "SNP calling " << files[i] << " ...\n";
                ind2species.push_back(name2id[meta.mappedname(indNames[i])]);
                SeqParser seq(files[i]);
                while (seq.nextSeq()){
                    string seqs = seq.getSeq(intqcs), seqn = seq.getSeq(intqcn);
                    snp.add(seqs, seqn);
                }

                for (const char c: snp.get()) tripInit.seq.append(c);
                snp.postprocess();
            }
        }
        SNP<K>::staticClear();
        */

       delete pSNP;
    }
};

int main(int argc, char** argv){
    ARG.setProgramName("waster", "Without-Alignment/Assembly Species Tree EstimatoR †");
    ARG.addIntArg(0, "sampled", 16, "Maximum number of sampled species for generating frequent patterns");
    ARG.addIntArg(0, "mode", 1, "1 (default): run the whole inferece, 2: only generate frequent patterns, 3: only generate SNPs, 4: start from SNPs");
    ARG.addStringArg(0, "continue", "", "Continue from provided frequent patterns");
    ARG.addIntArg(0, "qcs", 30, "Quality control threshold for the SNP bases (between 0-93, 30 by default)");
    ARG.addIntArg(0, "qcn", 20, "Quality control threshold for non-SNP bases (between 0-93, 20 by default)");
    ARG.addIntArg(0, "pattern", 500000000, "The number of most frequent k-mer patterns selected to generate consensus profiles");
    ARG.addIntArg(0, "consensus", 25000000, "The number of selected cosensus profiles to generate SNP table");

    Workflow WF(argc, argv);
    LOG << "#Base: " << WF.meta.tripInit.seq.len() << endl;
    auto res = WF.meta.run();
    LOG << "Normalized score: " << (double) res.first / 4 << endl;
    
    return 0;
}

#define DRIVER_VERSION "0"

#include<iostream>
#include<fstream>
#include<unordered_map>
#include<unordered_set>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<algorithm>
#include<random>
#include<thread>
#include<mutex>

#define ROOTING
#define NAME_MAPPING
//#define BLOCK_BOOTSTRAP

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

unsigned char CHAR2BITS[256] = {}, CHAR2MASK[256] = {}, RC[256] = {};

const string STRING_REVERSE_COMPLEMENT(const string &seq){
    string res;
    for (int i = seq.size() - 1; i >= 0; i--){
        res += RC[seq[i]];
    }
    return res;
}

template<size_t K> string DISPLAY(size_t bits){
    string res0, res1;
    for (int i = K * 2 - 1; i >= K; i--){
        size_t b = (bits >> (2 * i)) & 3;
        if (b == 0) res0 += 'A';
        if (b == 1) res0 += 'C';
        if (b == 2) res0 += 'G';
        if (b == 3) res0 += 'T';
    }
    for (int i = K - 1; i >= 0; i--){
        size_t b = (bits >> (2 * i)) & 3;
        if (b == 0) res1 += 'A';
        if (b == 1) res1 += 'C';
        if (b == 2) res1 += 'G';
        if (b == 3) res1 += 'T';
    }
    return res0 + "-" + res1 + "/" + STRING_REVERSE_COMPLEMENT(res1) + "-" + STRING_REVERSE_COMPLEMENT(res0);
}

template<size_t K> class KmerTable{
    constexpr static const size_t MASK = (1LL << (K * 2)) - 1;
    constexpr static const size_t LEN = 1LL << (K * 4);

    uint8_t* table;
public:

    static void shift(size_t &bits, size_t &mask, const unsigned char c){
        bits = ((bits << 2) | CHAR2BITS[c]) & MASK;
        mask = ((mask << 2) | CHAR2MASK[c]) & MASK;
    }

    static size_t REVERSE_COMPLEMENT(size_t bits){
        size_t res = 0;
        for (size_t i = 0; i < 2 * K; i++){
            res = (res << 2) | (3 ^ (bits & 3));
            bits >>= 2;
        }
        return res;
    }

    KmerTable(): table(new uint8_t[LEN]{}){}

    ~KmerTable(){
        delete table;
    }

    void add(const string &seq, bool bothStrands = true){
        size_t b1 = 0, b2 = 0, b3 = 0, b4 = 0;
        size_t m1 = 0, m2 = 0, m3 = 0, m4 = 0;
        for (size_t i = 0; i + 3 * K + 1 < seq.size(); i++){
            shift(b1, m1, seq[i + 0]);
            shift(b2, m2, seq[i + K]);
            shift(b3, m3, seq[i + 2 * K + 1]);
            shift(b4, m4, seq[i + 3 * K + 1]);
            if (i >= K - 1 && CHAR2BITS[seq[i + K + 1]] < 2 && CHAR2MASK[seq[i + K + 1]] == 0){
                if (m1 == 0 && m4 == 0) table[(b1 << (K * 2)) | b4] |= 1;
                if (m2 == 0 && m3 == 0) table[(b2 << (K * 2)) | b3] |= 1;
            }
        }
        if (bothStrands) add(STRING_REVERSE_COMPLEMENT(seq), false);
    }

    void postprocess(){
        for (size_t i = 0; i < LEN; i++){
            if (table[i] & 1){
                size_t j = REVERSE_COMPLEMENT(i);
                if (table[j] & 1){
                    if ((table[i] & 6) == 6) table[i] -= 1;
                    else table[i] += 1;
                    if (i == j) continue;
                    if ((table[j] & 6) == 6) table[j] -= 1;
                    else table[j] += 1;
                }
                else {
                    if ((table[i] & 248) == 248) table[i] -= 1;
                    else table[i] += 7;
                }
            }
        }
    }

    vector<size_t> frequentPatterns(size_t n) const{
        vector<size_t> freq[64];
        size_t cnt = 0, threshold = 2;
        for (size_t i = 0; i < LEN; i++){
            if (table[i] == 0) continue;
            if (table[i] & 6) continue;
            size_t j = REVERSE_COMPLEMENT(i);
            if (j <= i || table[j] == 0) continue;
            size_t level = (table[i] >> 3) + (table[j] >> 3);
            if (level < threshold) continue;
            freq[level].push_back(i);
            cnt++;
            while (cnt >= n) {
                cnt -= freq[threshold].size();
                threshold++;
            }
        }
        vector<size_t> res;
        for (size_t level = 63; level > 0; level--){
            for (size_t i: freq[level]){
                if (res.size() >= n) break;
                res.push_back(i);
                //cerr << DISPLAY<K>(i) << endl;
            }
        }
        return res;
    }
};

template<size_t K> class SNP{
    unordered_map<size_t, size_t> pattern2pos;
    string snps;

public:
    SNP(const vector<size_t> freqPatterns): pattern2pos(freqPatterns.size()), snps(freqPatterns.size(), (char) 0){
        for (int i = 0; i < freqPatterns.size(); i++){
            pattern2pos[freqPatterns[i]] = i;
        }
    }

    void add(size_t pattern, char c){
        if (pattern2pos.count(pattern) == 0) return;
        size_t pos = pattern2pos.at(pattern);
        if (snps[pos] == 0) snps[pos] = c;
        else if (snps[pos] != c) snps[pos] = 'N';
    }

    void add(const string &seq, bool bothStrands = true){
        size_t b1 = 0, b2 = 0, b3 = 0, b4 = 0;
        size_t m1 = 0, m2 = 0, m3 = 0, m4 = 0;
        for (size_t i = 0; i + 3 * K + 1 < seq.size(); i++){
            KmerTable<K>::shift(b1, m1, seq[i + 0]);
            KmerTable<K>::shift(b2, m2, seq[i + K]);
            KmerTable<K>::shift(b3, m3, seq[i + 2 * K + 1]);
            KmerTable<K>::shift(b4, m4, seq[i + 3 * K + 1]);
            if (i >= K - 1 && CHAR2MASK[seq[i + K + 1]] == 0){
                if (m1 == 0 && m4 == 0) add((b1 << (K * 2)) | b4, seq[i + K + 1]);
                if (m2 == 0 && m3 == 0) add((b2 << (K * 2)) | b3, seq[i + K + 1]);
            }
        }
        if (bothStrands) add(STRING_REVERSE_COMPLEMENT(seq), false);
    }

    void reset(){
        for (char& c: snps) c = 0;
    }

    const string& get(){
        for (char& c: snps){
            if (c == 0) c = 'N';
        }
        return snps;
    }
};

struct Workflow {
    MetaAlgorithm meta;
    TripartitionInitializer& tripInit = meta.tripInit;

    vector<string>& names = meta.names;
    unordered_map<string, int>& name2id = meta.name2id;

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
        gene.pi[0] = 1 - GCcontent;
        gene.pi[1] = GCcontent;
        gene.pi[2] = GCcontent;
        gene.pi[3] = 1 - GCcontent;
    }

    void formatGene(size_t pos, size_t nSite, size_t offset) {
        size_t nInd = ind2species.size(), nSpecies = names.size(), nKernal = nSite, nRep = 0;
        TripartitionInitializer::Gene::Initializer gene(nInd, nSpecies, nSite, nKernal, nRep);
        buildGeneSeq(gene, pos, nSite, offset);
        tripInit.genes.emplace_back(gene);
    }

    void readFasta(){
        ifstream fin(ARG.getStringArg("input"));
        size_t freqAT = 0, freqCG = 0, npos = 0;
        string line;
        while (getline(fin, line)){
            if (line[0] == '>'){
                string name = line.substr(1);
                LOG << "Processing " << name << " ...\n";
                addName(meta.mappedname(name));
                ind2species.push_back(name2id[meta.mappedname(name)]);
                continue;
            }
            npos += line.size();
            for (const char c: line) {
                tripInit.seq.append(c);
                switch (c) {
                    case 'A': case 'a': freqAT++; break;
                    case 'C': case 'c': freqCG++; break;
                    case 'G': case 'g': freqCG++; break;
                    case 'T': case 't': freqAT++; break;
                    case 'U': case 'u': freqAT++; break;
                }
            }
        }

        GCcontent = ((double) freqCG) / (freqAT + freqCG);
        nSNP = npos / ind2species.size();
    }

    Workflow(int argc, char** argv){
        //string mappingFile;
        meta.initialize(argc, argv);
        if (ARG.getStringArg("root") != "") addName(ARG.getStringArg("root"));
        tripInit.nThreads = meta.nThreads;

        if (ARG.getIntArg("mode") <= 3){
            switch (ARG.getIntArg("kmer")){
                case 7: init<7>(); break;
                case 8: init<8>(); break;
                case 9: init<9>(); break;
                default: cerr << "Bad k-mer size!\n"; exit(0);
            }
        }
        else readFasta();

        size_t nChunk = meta.nThreads;
        size_t pos = 0, len = nSNP;
        for (int i = 0; i < nChunk; i++) {
            size_t s = i * len / nChunk, t = (i + 1) * len / nChunk;
            formatGene(pos + s, t - s, len);
        }

        tripInit.nSpecies = names.size();
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
        }
        vector<size_t> freqPatterns;
        if (ARG.getStringArg("continue") == ""){
            unordered_set<size_t> selected;
            KmerTable<K> table;
            size_t nSample = ARG.getIntArg("sampled");
            for (size_t i = 0; i < files.size() && i < nSample; i++){
                size_t cur = rand() % files.size();
                while (selected.count(cur)) cur = rand() % files.size();
                selected.insert(cur);
                LOG << "Species " << indNames[cur] << " is selected to count the most frequent patterns.\n";
                ifstream fin(files[cur]);
                string seq, line;
                while (getline(fin, line)){
                    if (line[0] == '>'){
                        if (seq != "") table.add(seq);
                        seq = "";
                    }
                    else{
                        seq += line;
                    }
                }
                if (seq != "") table.add(seq);
                table.postprocess();
            }
            freqPatterns = table.frequentPatterns((4LL << (4 * K)) / files.size());

            if (ARG.getIntArg("mode") == 2){
                if (ARG.getStringArg("output") == "<standard output>"){
                    for (size_t i: freqPatterns) cout << i << endl;
                }
                else {
                    ofstream fout(ARG.getStringArg("output"));
                    for (size_t i: freqPatterns) fout << i << endl;
                }
                exit(0);
            }
        }
        else {
            ifstream fin(ARG.getStringArg("continue"));
            size_t i;
            while (fin >> i) freqPatterns.push_back(i);
        }
        
        nSNP = freqPatterns.size();
        LOG << nSNP << " SNPs are selected.\n";
        
        if (ARG.getIntArg("mode") == 3){
            if (ARG.getStringArg("output") != "<standard output>") ofstream fin(ARG.getStringArg("output"));
        }

        SNP<K> snp(freqPatterns);
        freqPatterns.clear();
        size_t freqAT = 0, freqCG = 0;
        for (size_t i = 0; i < files.size(); i++){
            LOG << "Processing " << files[i] << " ...\n";
            ind2species.push_back(name2id[meta.mappedname(indNames[i])]);
            ifstream fin(files[i]);
            string seq, line;
            while (getline(fin, line)){
                if (line[0] == '>'){
                    if (seq != "") snp.add(seq);
                    seq = "";
                }
                else{
                    seq += line;
                }
            }
            if (seq != "") snp.add(seq);
            
            if (ARG.getIntArg("mode") == 3){
                if (ARG.getStringArg("output") == "<standard output>"){
                    cout << ">" << indNames[i] << endl << snp.get() << endl;
                }
                else {
                    ofstream fout(ARG.getStringArg("output"), ios_base::app);
                    fout << ">" << indNames[i] << endl << snp.get() << endl;
                }
            }

            for (const char c: snp.get()) {
                tripInit.seq.append(c);
                switch (c) {
                    case 'A': case 'a': freqAT++; break;
                    case 'C': case 'c': freqCG++; break;
                    case 'G': case 'g': freqCG++; break;
                    case 'T': case 't': freqAT++; break;
                    case 'U': case 'u': freqAT++; break;
                }
            }
            snp.reset();
        }

        if (ARG.getIntArg("mode") == 3) exit(0);
        GCcontent = ((double) freqCG) / (freqAT + freqCG);
    }
};

int main(int argc, char** argv){
    CHAR2BITS['a'] = 0; CHAR2BITS['c'] = 1; CHAR2BITS['g'] = 2; CHAR2BITS['t'] = 3; CHAR2BITS['u'] = 3;
    CHAR2BITS['A'] = 0; CHAR2BITS['C'] = 1; CHAR2BITS['G'] = 2; CHAR2BITS['T'] = 3; CHAR2BITS['U'] = 3;

    for (size_t i = 0; i < 256; i++) CHAR2MASK[i] = 1;
    CHAR2MASK['a'] = 0; CHAR2MASK['c'] = 0; CHAR2MASK['g'] = 0; CHAR2MASK['t'] = 0; CHAR2MASK['u'] = 0;
    CHAR2MASK['A'] = 0; CHAR2MASK['C'] = 0; CHAR2MASK['G'] = 0; CHAR2MASK['T'] = 0; CHAR2MASK['U'] = 0;

    RC['a'] = 't'; RC['c'] = 'g'; RC['g'] = 'c'; RC['t'] = 'a'; RC['u'] = 'a'; 
    RC['A'] = 'T'; RC['C'] = 'G'; RC['G'] = 'C'; RC['T'] = 'A'; RC['U'] = 'A'; 

    ARG.setProgramName("waster-site", "Without-Alignment/Assembly Species Tree EstimatoR † (site)");
    ARG.addIntArg('k', "kmer", 8, "k-mer size; 7: require >256 MB memory, 8 (default): >4 GB memory, 9: >64 GB memory", true);
    ARG.addIntArg(0, "sampled", 64, "Maximum number of sampled species for generating frequent patterns");
    ARG.addIntArg(0, "mode", 1, "1 (default): run the whole inferece, 2: only generate frequent patterns, 3: only generate SNPs, 4: start from SNPs");
    ARG.addStringArg(0, "continue", "", "Continue from provided frequent patterns");

    Workflow WF(argc, argv);
    LOG << "#Base: " << WF.meta.tripInit.seq.len() << endl;
    auto res = WF.meta.run();
    LOG << "Score: " << (double) res.first << endl;
    
    return 0;
}

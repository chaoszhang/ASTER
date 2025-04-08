#define DRIVER_VERSION "2"

/* CHANGE LOG
 * 2: Final score normalization
 * 1: Remove non-effective sites
 */

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

static unsigned char CHAR2BITS(char c){
    if (c == 'A') return 0;
    if (c == 'C') return 1;
    if (c == 'G') return 2;
    if (c == 'T') return 3;
    return 0;
}

static unsigned char CHAR2MASK(char c){
    if (c == 'N') return 1;
    return 0;
}

template<size_t K> class KmerTable{
    constexpr static const size_t MASK = (1LL << (K * 2)) - 1;

    uint8_t* table;

public:
    constexpr static const size_t LEN = (1LL << (K * 4 - 1)) - (1LL << (K * 2 - 1));

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

    KmerTable(): table(new uint8_t[LEN]{}){}

    ~KmerTable(){
        delete table;
    }

    void add(const string &seqs, const string &seqn, bool bothStrands = true){
        size_t b1 = 0, b2 = 0, b3 = 0, b4 = 0;
        size_t m1 = 0, m2 = 0, m3 = 0, m4 = 0;
        for (size_t i = 0; i + 3 * K + 1 < seqn.size(); i++){
            shift(b1, m1, seqn[i + 0]);
            shift(b2, m2, seqn[i + K]);
            rshift(b3, m3, seqn[i + 2 * K + 1]);
            rshift(b4, m4, seqn[i + 3 * K + 1]);
            char c = seqs[i + K + 1];
            if (i >= K - 1 && CHAR2MASK(c) == 0){
                if (b1 > b4 && m1 == 0 && m4 == 0){
                    size_t p = b1 * (b1 - 1) / 2 + b4;
                    if (CHAR2BITS(c) & 1) table[p] |= 2;
                    else table[p] |= 1;
                }
                if (b2 > b3 && m2 == 0 && m3 == 0){
                    size_t p = b2 * (b2 - 1) / 2 + b3;
                    if (CHAR2BITS(c) & 1) table[p] |= 2;
                    else table[p] |= 1;
                }
            }
        }
        if (bothStrands) add(SeqParser::STRING_REVERSE_COMPLEMENT(seqs), SeqParser::STRING_REVERSE_COMPLEMENT(seqn), false);
    }

    void postprocess(){
        fillcnt = 0;
        for (size_t i = 0; i < LEN; i++){
            if (table[i] == 0) continue; 
            fillcnt++;
            if ((table[i] & 3) == 0) continue;
            if ((table[i] & 3) == 3) continue;
            table[i] &= 252;
            if (table[i] < 252) table[i] += 4;
        }
    }

    double fillProportion() const{
        return fillcnt / (double) LEN;
    }

    vector<size_t> frequentPatterns(size_t n, bool optimize4noDownsample = false) const{
        vector<size_t> freq[64];
        size_t cnt = 0, threshold = (optimize4noDownsample) ? 4 : 2;
        for (size_t i = 0; i < LEN; i++){
            if (table[i] == 0) continue;
            if (table[i] & 3) continue;
            size_t level = table[i] / 4;
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
            }
        }
        return res;
    }
};

template<size_t K> class SNP{
    constexpr static const size_t SHIFT = 8;
    constexpr static const size_t M = (1LL << (K * 4 - 1 - SHIFT)) - (1LL << (K * 2 - 1 - SHIFT));
    string snps;

public:
    static vector<size_t> sortedPattern, pattern2pos;
    static vector<int8_t> cntAG, cntCT;

    static void staticInit(){
        sort(sortedPattern.begin(), sortedPattern.end());
        cntAG.resize(sortedPattern.size());
        cntCT.resize(sortedPattern.size());
        pattern2pos.reserve(M + 1);
        size_t n = sortedPattern.size();
        for (size_t i = 0, j = 0; i < M + 1; i++){
            while (j < n && sortedPattern[j] < (i << SHIFT)) j++;
            pattern2pos.push_back(j);
        }
    }

    static void staticFilter(){
        vector<size_t> filteredPattern;
        for (size_t i = 0; i < sortedPattern.size(); i++) {
            #ifdef CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH
            filteredPattern.push_back(sortedPattern[i]);
            #else
            if (cntAG[i] >= 2 && cntCT[i] >= 2) filteredPattern.push_back(sortedPattern[i]);
            #endif
        }
        sortedPattern = filteredPattern;
        pattern2pos.clear();
        cntAG.clear();
        cntCT.clear();
        staticInit();
    }

    static void staticClear(){
        sortedPattern.clear();
        sortedPattern.shrink_to_fit();
        pattern2pos.clear();
        pattern2pos.shrink_to_fit();
    }

    static size_t locate(size_t pattern){
        size_t L = pattern2pos[pattern >> SHIFT], R = pattern2pos[(pattern >> SHIFT) + 1];
        size_t i = L;
        while (i < R && sortedPattern[i] < pattern) i++;
        if (sortedPattern[i] == pattern) return i;
        return -1;
    }

    SNP(): snps(sortedPattern.size(), (char) 0){}

    void add(size_t pattern, char c){
        size_t pos = locate(pattern);
        if (pos == static_cast<size_t>(-1)) return;
        if (snps[pos] == 0) snps[pos] = c;
        else if (snps[pos] != c) snps[pos] = 'N';
    }

    void add(const string &seqs, const string &seqn, bool bothStrands = true){
        size_t b1 = 0, b2 = 0, b3 = 0, b4 = 0;
        size_t m1 = 0, m2 = 0, m3 = 0, m4 = 0;
        for (size_t i = 0; i + 3 * K + 1 < seqs.size(); i++){
            KmerTable<K>::shift(b1, m1, seqn[i + 0]);
            KmerTable<K>::shift(b2, m2, seqn[i + K]);
            KmerTable<K>::rshift(b3, m3, seqn[i + 2 * K + 1]);
            KmerTable<K>::rshift(b4, m4, seqn[i + 3 * K + 1]);
            if (i >= K - 1 && CHAR2MASK(seqs[i + K + 1]) == 0){
                if (b1 > b4 && m1 == 0 && m4 == 0) add(b1 * (b1 - 1) / 2 + b4, seqs[i + K + 1]);
                if (b2 > b3 && m2 == 0 && m3 == 0) add(b2 * (b2 - 1) / 2 + b3, seqs[i + K + 1]);
            }
        }
        if (bothStrands) add(SeqParser::STRING_REVERSE_COMPLEMENT(seqs), SeqParser::STRING_REVERSE_COMPLEMENT(seqn), false);
    }

    void postprocess(){
        for (size_t i = 0; i < sortedPattern.size(); i++) {
            if ((snps[i] == 'A' || snps[i] == 'G') && cntAG[i] < 2) cntAG[i]++;
            if ((snps[i] == 'C' || snps[i] == 'T') && cntCT[i] < 2) cntCT[i]++; 
            snps[i] = 0;
        } 
    }

    const string& get(){
        for (char& c: snps){
            if (c == 0) c = 'N';
        }
        return snps;
    }
};
template <size_t K> vector<size_t> SNP<K>::sortedPattern;
template <size_t K> vector<size_t> SNP<K>::pattern2pos;
template <size_t K> vector<int8_t> SNP<K>::cntAG;
template <size_t K> vector<int8_t> SNP<K>::cntCT;

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
            switch (ARG.getIntArg("kmer")){
                case 7: init<7>(); break;
                case 8: init<8>(); break;
                case 9: init<9>(); break;
                case 10: init<10>(); break;
                default: cerr << "Bad k-mer size!\n"; exit(0);
            }
        }
        else readAlignment();

        size_t nChunk = meta.nThreads;
        size_t pos = 0, len = nSNP;
        for (size_t i = 0; i < nChunk; i++) {
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
            if (files.size() == 0){
                cerr << "Error: Input list file '" << ARG.getStringArg("input") <<"' empty! Typo?\n";
                exit(0);
            }
        }
        vector<size_t> &freqPatterns = SNP<K>::sortedPattern;
        if (ARG.getStringArg("continue") == ""){
            unordered_set<size_t> selected;
            KmerTable<K> table;
            size_t nSample = ARG.getIntArg("sampled");
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
                LOG << "Hash table " << (int) (table.fillProportion() * 100) << "% filled." << endl;
                if (i >= 3 && table.fillProportion() > 0.33){
                    LOG << "Early termination of k-mer search due to high fill rate after processing " << i + 1 << " samples. You may consider increasing K for a larger Hash table if this number is too low." << endl;
                    break;
                }
            }
            size_t patternCnt = KmerTable<K>::LEN * 2 / (files.size() + 128);
            freqPatterns = table.frequentPatterns(patternCnt, selected.size() == files.size());

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
        
        LOG << freqPatterns.size() << " SNPs are selected.\n";
        
        if (ARG.getIntArg("mode") == 3){
            if (ARG.getStringArg("output") != "<standard output>") ofstream fin(ARG.getStringArg("output"));
        }
        
        SNP<K>::staticInit();
        {
            SNP<K> snp;
            size_t freqAT = 0, freqCG = 0;
            for (size_t i = 0; i < files.size(); i++){
                LOG << "Processing " << files[i] << " ...\n";
                SeqParser seq(files[i]);
                while (seq.nextSeq()){
                    string seqs = seq.getSeq(intqcs), seqn = seq.getSeq(intqcn);
                    snp.add(seqs, seqn);
                }
                if (ARG.getIntArg("mode") == 3){
                    if (ARG.getStringArg("output") == "<standard output>"){
                        cout << ">" << indNames[i] << endl << snp.get() << endl;
                    }
                    else {
                        ofstream fout(ARG.getStringArg("output"), ios_base::app);
                        fout << ">" << indNames[i] << endl << snp.get() << endl;
                    }
                }

                string s = snp.get();
                freqAT += countAT(s);
                freqCG += countCG(s);
                snp.postprocess();
            }
            if (ARG.getIntArg("mode") == 3) exit(0);
            GCcontent = ((double) freqCG) / (freqAT + freqCG);
        }

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
    }
};

int main(int argc, char** argv){
    ARG.setProgramName("waster-site", "Without-Alignment/Assembly Species Tree EstimatoR † (site)\n***Running WASTER requires 64GB memory by default!***");
    ARG.addIntArg('k', "kmer", 9, "k-mer size; 7: require >128 MB memory, 8: >2 GB memory, 9 (default): >32 GB memory, 10: >512 GB memory", true);
    ARG.addIntArg(0, "sampled", 64, "Maximum number of sampled species for generating frequent patterns");
    ARG.addIntArg(0, "mode", 1, "1 (default): run the whole inferece, 2: only generate frequent patterns, 3: only generate SNPs, 4: start from SNPs");
    ARG.addStringArg(0, "continue", "", "Continue from provided frequent patterns");
    ARG.addIntArg(0, "qcs", 30, "Quality control threshold for the SNP base (between 0-93, 30 by default)");
    ARG.addIntArg(0, "qcn", 20, "Quality control threshold for non-SNP bases (between 0-93, 20 by default)");

    Workflow WF(argc, argv);
    LOG << "#Base: " << WF.meta.tripInit.seq.len() << endl;
    auto res = WF.meta.run();
    LOG << "Normalized score: " << (double) res.first / 4 << endl;
    
    return 0;
}

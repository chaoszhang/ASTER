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
unsigned char QUALITY2ASCII[95] = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

const string STRING_REVERSE_COMPLEMENT(const string &seq){
    string res;
    for (int i = seq.size() - 1; i >= 0; i--){
        res += RC[seq[i]];
    }
    return res;
}

template<size_t K> class KmerTable{
    constexpr static const size_t MASK = (1LL << (K * 2)) - 1;

    uint8_t* table;

public:
    constexpr static const size_t LEN = (1LL << (K * 4 - 1)) - (1LL << (K * 2 - 1));

    size_t fillcnt = 0;

    static void shift(size_t &bits, size_t &mask, const unsigned char c){
        bits = ((bits << 2) | CHAR2BITS[c]) & MASK;
        mask = ((mask << 2) | CHAR2MASK[c]) & MASK;
    }

    static void rshift(size_t &bits, size_t &mask, const unsigned char c){
        bits = (bits >> 2) | ((3 ^ CHAR2BITS[c]) << (K * 2 - 2));
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

    void add(const string &seqs, const string &seqn, bool bothStrands = true){
        size_t b1 = 0, b2 = 0, b3 = 0, b4 = 0;
        size_t m1 = 0, m2 = 0, m3 = 0, m4 = 0;
        for (size_t i = 0; i + 3 * K + 1 < seqn.size(); i++){
            shift(b1, m1, seqn[i + 0]);
            shift(b2, m2, seqn[i + K]);
            rshift(b3, m3, seqn[i + 2 * K + 1]);
            rshift(b4, m4, seqn[i + 3 * K + 1]);
            if (i >= K - 1 && CHAR2MASK[seqs[i + K + 1]] == 0){
                if (b1 > b4 && m1 == 0 && m4 == 0){
                    size_t p = b1 * (b1 - 1) / 2 + b4;
                    if (CHAR2BITS[seqs[i + K + 1]] & 1) table[p] |= 2;
                    else table[p] |= 1;
                }
                if (b2 > b3 && m2 == 0 && m3 == 0){
                    size_t p = b2 * (b2 - 1) / 2 + b3;
                    if (CHAR2BITS[seqs[i + K + 1]] & 1) table[p] |= 2;
                    else table[p] |= 1;
                }
            }
        }
        if (bothStrands) add(STRING_REVERSE_COMPLEMENT(seqs), STRING_REVERSE_COMPLEMENT(seqn), false);
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
    constexpr static const size_t SHIFT = 6;
    constexpr static const size_t M = (1LL << (K * 4 - 1 - SHIFT)) - (1LL << (K * 2 - 1 - SHIFT));
    string snps;

public:
    static vector<size_t> sortedPattern, pattern2pos;

    static void staticInit(){
        sort(sortedPattern.begin(), sortedPattern.end());
        pattern2pos.reserve(M + 1);
        size_t n = sortedPattern.size();
        for (size_t i = 0, j = 0; i < M + 1; i++){
            while (j < n && sortedPattern[j] < (i << SHIFT)) j++;
            pattern2pos.push_back(j);
        }
    }

    static void staticClear(){
        sortedPattern.clear();
        pattern2pos.clear();
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
        if (pos == -1) return;
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
            if (i >= K - 1 && CHAR2MASK[seqs[i + K + 1]] == 0){
                if (b1 > b4 && m1 == 0 && m4 == 0) add(b1 * (b1 - 1) / 2 + b4, seqs[i + K + 1]);
                if (b2 > b3 && m2 == 0 && m3 == 0) add(b2 * (b2 - 1) / 2 + b3, seqs[i + K + 1]);
            }
        }
        if (bothStrands) add(STRING_REVERSE_COMPLEMENT(seqs), STRING_REVERSE_COMPLEMENT(seqn), false);
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
template <size_t K> vector<size_t> SNP<K>::sortedPattern;
template <size_t K> vector<size_t> SNP<K>::pattern2pos;

struct Workflow {
    MetaAlgorithm meta;
    TripartitionInitializer& tripInit = meta.tripInit;

    vector<string>& names = meta.names;
    unordered_map<string, int>& name2id = meta.name2id;
    
    unsigned char QCS = 0, QCN = 0;
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

    static string getFormat(const string fileName){
        ifstream fin(fileName);
        string s;
        getline(fin, s);
        if (s.size() == 0) return "";
        if (s[0] == '@') return "fastq";
        if (s[0] == '>') return "fasta";
        return "";
    }

    bool getFastaSeq(ifstream& fin, string& seq){
        if (eofFlag){
            eofFlag = false;
            return false;
        }
        seq = "";
        eofFlag = true;
        string line;
        while(getline(fin, line)){
            if (line[0] == '>'){
                if (seq == "") continue;
                else {
                    eofFlag = false;
                    break;
                }
            }
            seq += line;
        }
        return true;
    }

    bool getFastqSeq(ifstream& fin, string& seqs, string& seqn){
        string line;
        if (getline(fin, line)){
            if (line == "") getline(fin, line);
            getline(fin, seqs);
            seqn = seqs;
            getline(fin, line);
            getline(fin, line);
            for (int i = 0; i < seqs.size(); i++){
                if (line[i] < QCS) seqs[i] = 'N';
                if (line[i] < QCN) seqn[i] = 'N';
            }

            return true;
        }
        else return false;
    }

    void readAlignment(){
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
        LOG << "Make sure you have run 'waster-site -h', read about '-k' command, and ensured you have enough memory to proceed!\n";
        if (ARG.getStringArg("root") != "") addName(ARG.getStringArg("root"));
        tripInit.nThreads = meta.nThreads;

        if (ARG.getIntArg("mode") <= 3){
            int intqcs = ARG.getIntArg("qcs");
            if (intqcs >= 94){
                cerr << "Bad quality control threshold for SNPs!\n";
                exit(0);
            }
            QCS = QUALITY2ASCII[intqcs];
            LOG << "Quality control: Masking all SNP bases with quality lower than '" << QCS << "' for FASTQ inputs.\n";
            int intqcn = ARG.getIntArg("qcn");
            if (intqcn >= 94){
                cerr << "Bad quality control threshold for non-SNPs!\n";
                exit(0);
            }
            QCN = QUALITY2ASCII[intqcn];
            LOG << "Quality control: Masking all non-SNP bases with quality lower than '" << QCN << "' for FASTQ inputs.\n";
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
                string format = getFormat(files[cur]);
                ifstream fin(files[cur]);
                if (format == "fasta"){
                    string seq;
                    while (getFastaSeq(fin, seq)) table.add(seq, seq);
                }
                else if (format == "fastq"){
                    string seqs, seqn;
                    while (getFastqSeq(fin, seqs, seqn)) table.add(seqs, seqn);
                }
                else {
                    cerr << "File " << files[cur] << " bad format!\n";
                    exit(0);
                }
                table.postprocess();
                LOG << "Hash table " << (int) (table.fillProportion() * 100) << "% filled." << endl;
                if (i >= 3 && table.fillProportion() > 0.33){
                    LOG << "Early termination of k-mer search due to high fill rate after processing" << i + 1 << " samples. You may consider increasing K for a larger Hash table if this number is too low." << endl;
                    break;
                }
            }
            size_t patternCnt = min(KmerTable<K>::LEN * 4 / files.size(), KmerTable<K>::LEN / 8);
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
        
        nSNP = freqPatterns.size();
        LOG << nSNP << " SNPs are selected.\n";
        
        if (ARG.getIntArg("mode") == 3){
            if (ARG.getStringArg("output") != "<standard output>") ofstream fin(ARG.getStringArg("output"));
        }

        SNP<K>::staticInit();
        SNP<K> snp;
        //freqPatterns.clear();
        size_t freqAT = 0, freqCG = 0;
        for (size_t i = 0; i < files.size(); i++){
            LOG << "Processing " << files[i] << " ...\n";
            ind2species.push_back(name2id[meta.mappedname(indNames[i])]);
            string format = getFormat(files[i]);
            ifstream fin(files[i]);
            if (format == "fasta"){
                string seq;
                while (getFastaSeq(fin, seq)) snp.add(seq, seq);
            }
            else if (format == "fastq"){
                string seqs, seqn;
                while (getFastqSeq(fin, seqs, seqn)) snp.add(seqs, seqn);
            }
            else {
                cerr << "File " << files[i] << " bad format!\n";
                exit(0);
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
        SNP<K>::staticClear();
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
    ARG.addIntArg('k', "kmer", 9, "k-mer size; 7: require >128 MB memory, 8: >2 GB memory, 9 (default): >32 GB memory, 10: >512 GB memory", true);
    ARG.addIntArg(0, "sampled", 64, "Maximum number of sampled species for generating frequent patterns");
    ARG.addIntArg(0, "mode", 1, "1 (default): run the whole inferece, 2: only generate frequent patterns, 3: only generate SNPs, 4: start from SNPs");
    ARG.addStringArg(0, "continue", "", "Continue from provided frequent patterns");
    ARG.addIntArg(0, "qcs", 30, "Quality control threshold for the SNP base (between 0-93, 30 by default)");
    ARG.addIntArg(0, "qcn", 20, "Quality control threshold for non-SNP bases (between 0-93, 20 by default)");

    Workflow WF(argc, argv);
    LOG << "#Base: " << WF.meta.tripInit.seq.len() << endl;
    auto res = WF.meta.run();
    LOG << "Score: " << (double) res.first << endl;
    
    return 0;
}

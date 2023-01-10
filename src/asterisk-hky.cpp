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

int GENE_ID = 0;

struct Workflow {
    MetaAlgorithm meta;
    TripartitionInitializer& tripInit = meta.tripInit;

    vector<string>& names = meta.names;
    unordered_map<string, int>& name2id = meta.name2id;

    string guide;
    int diskcover;
    bool useDCM = false;

    void addName(const string& name) {
        if (name2id.count(name) == 0) {
            name2id[name] = names.size();
            names.push_back(name);
        }
    }

    static int infoCount(const array<int, 4> &cnt){
        int result = 0;
        for (int j = 0; j < 4; j++) {
            if (cnt[j] >= 2) result++;
        }
        return result;
    }

    static array<int, 4> add(const array<int, 4> &a, const array<int, 4> &b){
        array<int, 4> result;
        for (int j = 0; j < 4; j++) {
            result[j] = a[j] + b[j];
        }
        return result;
    }

    static int sum(const array<int, 4> &cnt){
        int result = 0;
        for (int j = 0; j < 4; j++) {
            result += cnt[j];
        }
        return result;
    }

    template<int THRESHOLD> array<int, 4> breakHierarchy(unordered_map<int, int> &result, 
            BinaryTree &hierarchy, int cur, size_t pos, size_t offset){
        if (hierarchy.taxon(cur) != -1){
            array<int, 4> curCnt;
            for (int j = 0; j < 4; j++){
                curCnt[j] = tripInit.seq.get(pos + hierarchy.taxon(cur) * offset, j);
            }
            return curCnt;
        }
        int lc = hierarchy.left(cur), rc = hierarchy.right(cur);
        array<int, 4> lv = breakHierarchy<THRESHOLD>(result, hierarchy, lc, pos, offset);
        array<int, 4> rv = breakHierarchy<THRESHOLD>(result, hierarchy, rc, pos, offset);
        array<int, 4> v = add(lv, rv);
        int cnt = result.size();
        if (infoCount(v) > THRESHOLD && hierarchy.breaking(cur)){
            if (!(infoCount(lv) > THRESHOLD && hierarchy.breaking(lc)) && sum(lv) >= 4 && infoCount(lv) > 1) result[lc] = cnt++;
            if (!(infoCount(rv) > THRESHOLD && hierarchy.breaking(rc)) && sum(rv) >= 4 && infoCount(rv) > 1) result[rc] = cnt++;
        }
        else {
            if (cur == hierarchy.root && sum(v) >= 4) result[cur] = cnt++;
        }
        return v;
    }

    void ind2kernal(vector<int> &result, const unordered_map<int, int> &kernals, BinaryTree &hierarchy, int cur, int kernal = -1){
        if (kernals.count(cur)) kernal = kernals.at(cur);
        if (hierarchy.taxon(cur) != -1){
            result[hierarchy.taxon(cur)] = kernal;
        }
        else {
            int lc = hierarchy.left(cur), rc = hierarchy.right(cur);
            ind2kernal(result, kernals, hierarchy, lc, kernal);
            ind2kernal(result, kernals, hierarchy, rc, kernal);
        }
    }

    void buildGeneSeq(TripartitionInitializer::Gene::Initializer &gene, const vector<int> &ind2species, size_t pos, size_t nSite, size_t offset){
        int nInd = ind2species.size();
        array<size_t, 4> cnt = {};
        for (int iInd = 0; iInd < nInd; iInd++) {
            size_t pSeq = pos + iInd * offset;
            for (int iSite = 0; iSite < nSite; iSite++) {
                tripInit.seq.add(pSeq + iSite, cnt);
            }
            gene.species2ind[ind2species[iInd]].push_back(iInd);
            gene.ind2seq[iInd] = pSeq;
        }
        double total = cnt[0] + cnt[1] + cnt[2] + cnt[3];
        for (int j = 0; j < 4; j++) {
            gene.pi[j] = cnt[j] / total;
        }
    }

    void formatGeneNaive(const vector<int> &ind2species, size_t pos, size_t nSite, size_t offset) {
        int nInd = ind2species.size(), nSpecies = names.size(), nKernal = nSite, nRep = 0;
        TripartitionInitializer::Gene::Initializer gene(nInd, nSpecies, nSite, nKernal, nRep);
        buildGeneSeq(gene, ind2species, pos, nSite, offset);
        for (int iKernal = 0; iKernal < nKernal; iKernal++){
            gene.weight[iKernal] = 1;
        }
        tripInit.genes.emplace_back(gene);
    }

    void formatGeneDCM(const vector<int> &ind2species, size_t pos, size_t nSite, size_t offset) {
        LOG << "Building local tree " << ++GENE_ID << " ...\n";
        Workflow WF(*this, ind2species, pos, nSite, offset);
        BinaryTree gtree = WF.getGenetree();
        vector<vector<vector<int> > > siteRepInd2kernal;
        vector<vector<int> > kernalOffset;
        int nInd = ind2species.size(), nSpecies = names.size(), nKernal = 0, nRep = diskcover;
        for (int iSite = 0; iSite < nSite; iSite++){
            siteRepInd2kernal.emplace_back();
            kernalOffset.emplace_back();
            for (int iRep = 0; iRep < nRep; iRep++){
                BinaryTree hierarchy = gtree.sample(); 
                unordered_map<int, int> kernals;
                breakHierarchy<3>(kernals, hierarchy, hierarchy.root, pos, offset);
                kernalOffset[iSite].push_back(nKernal);
                nKernal += kernals.size();
                siteRepInd2kernal[iSite].emplace_back(nInd);
                ind2kernal(siteRepInd2kernal[iSite][iRep], kernals, hierarchy, hierarchy.root, -1);
            }
        }
        TripartitionInitializer::Gene::Initializer gene(nInd, nSpecies, nSite, nKernal, nRep);
        buildGeneSeq(gene, ind2species, pos, nSite, offset);
        vector<int> kernalSize(nKernal);
        //ofstream test("test"); 
        for (int iInd = 0; iInd < nInd; iInd++) {
            for (int iSite = 0; iSite < nSite; iSite++) {
                for (int iRep = 0; iRep < nRep; iRep++){
                    int iKernal = siteRepInd2kernal[iSite][iRep][iInd];
                    if (iKernal != -1) iKernal += kernalOffset[iSite][iRep];
                    gene.setIndSiteRep2kernal(iInd, iSite, iRep, iKernal);
                    if (iKernal != -1) kernalSize[iKernal]++;
                    //test << iKernal << "\t";
                }
            }
            //test << endl;
        }
        for (int iKernal = 0; iKernal < nKernal; iKernal++){
            if (kernalSize[iKernal] < 3) gene.weight[iKernal] = 0;
            else gene.weight[iKernal] = 1;
        }
        tripInit.genes.emplace_back(gene);
    }

    void formatGene(const vector<int> &ind2species, size_t pos, size_t nSite, size_t offset) {
        if (useDCM) formatGeneDCM(ind2species, pos, nSite, offset);
        else formatGeneNaive(ind2species, pos, nSite, offset);
    }

    void readFasta(const string file) {
        vector<bool> keep;
        {
            vector<array<unsigned short, 4> > freq;
            ifstream fin(file);
            string name;
            fin >> name;
            while (name != "") {
                addName(name.substr(1));
                name = "";
                string seq, line;
                while (fin >> line) {
                    if (line[0] == '>') { name = line; break; }
                    seq += line;
                }
                if (freq.size() != seq.size()) {
                    if (freq.size() == 0) freq.resize(seq.size());
                    else { cerr << "File '" << file << "' is ill-formated."; exit(0); }
                }
                for (size_t i = 0; i < seq.size(); i++) {
                    switch (seq[i]) {
                        case 'A': case 'a': freq[i][0]++; break;
                        case 'C': case 'c': freq[i][1]++; break;
                        case 'G': case 'g': freq[i][2]++; break;
                        case 'T': case 't': freq[i][3]++; break;
                    }
                }
            }
            keep.resize(freq.size());
            for (size_t i = 0; i < keep.size(); i++) {
                int cnt = 0;
                for (int j = 0; j < 4; j++) {
                    if (freq[i][j] >= 2) cnt++;
                }
                keep[i] = (cnt >= 2);
            }
        }
        {
            ifstream fin(file);
            string name;
            fin >> name;
            vector<int> ind2species;
            size_t pos = tripInit.seq.len(), len = count(keep.begin(), keep.end(), true);
            while (name != "") {
                ind2species.push_back(name2id[name.substr(1)]);
                name = "";
                string line, seq;
                while (fin >> line) {
                    if (line[0] == '>') { name = line; break; }
                    seq += line;
                }
                for (int i = 0; i < keep.size(); i++){
                    if (keep[i]) tripInit.seq.append(seq[i]);
                }
            }
            int nChunk = (len + ARG.getIntArg("chunk") - 1) / ARG.getIntArg("chunk");
            for (int i = 0; i < nChunk; i++) {
                size_t s = i * len / nChunk, t = (i + 1) * len / nChunk;
                formatGene(ind2species, pos + s, t - s, len);
            }
        }
    }

    bool readPhylip(ifstream &fin, ifstream &fin2){
        size_t nTaxa, nSites;
        if (!(fin >> nTaxa)) return false;
        fin >> nSites;
        vector<bool> keep(nSites);
        {
            vector<array<unsigned short, 4> > freq(nSites);
            for (int i = 0; i < nTaxa; i++){
                string name, seq;
                fin >> name >> seq;
                if (seq.size() != nSites) { cerr << "The input is ill-formated."; exit(0); }
                addName(name);
                for (int j = 0; j < seq.size(); j++) {
                    switch (seq[j]) {
                        case 'A': case 'a': freq[j][0]++; break;
                        case 'C': case 'c': freq[j][1]++; break;
                        case 'G': case 'g': freq[j][2]++; break;
                        case 'T': case 't': freq[j][3]++; break;
                    }
                }
            }
            for (int i = 0; i < keep.size(); i++) {
                int cnt = 0;
                for (int j = 0; j < 4; j++) {
                    if (freq[i][j] >= 2) cnt++;
                }
                keep[i] = (cnt >= 2);
            }
        }
        {
            fin2 >> nTaxa >> nSites;
            size_t pos = tripInit.seq.len(), len = count(keep.begin(), keep.end(), true);
            vector<int> ind2species;
            for (int i = 0; i < nTaxa; i++){
                string name, seq;
                fin2 >> name >> seq;
                ind2species.push_back(name2id[name]);
                for (size_t j = 0; j < keep.size(); j++){
                    if (keep[j]) tripInit.seq.append(seq[j]);
                }
            }
            int nChunk = (len + ARG.getIntArg("chunk") - 1) / ARG.getIntArg("chunk");
            for (int i = 0; i < nChunk; i++) {
                size_t s = i * len / nChunk, t = (i + 1) * len / nChunk;
                formatGene(ind2species, pos + s, t - s, len);
            }
        }
        return true;
    }

    Workflow(int argc, char** argv){
        //string mappingFile;
        meta.initialize(argc, argv);
        init();
    }

    Workflow(int argc, char** argv, string guide, const vector<string> names, unordered_map<string, int> name2id):
            guide(guide), diskcover(ARG.getIntArg("diskcover")), useDCM(true) {
        //string mappingFile;
        meta.initialize(argc, argv);
        meta.names = names;
        meta.name2id = name2id;
        meta.nRounds = 0;
        meta.guideTree = guide;
        init();
    }

    void init(){
        tripInit.nThreads = meta.nThreads;
        
        if (ARG.getStringArg("format") == "fasta") {
            readFasta(ARG.getStringArg("input"));
        }
        else if (ARG.getStringArg("format") == "list") {
            ifstream fin(ARG.getStringArg("input"));
            string line;
            while (getline(fin, line)) {
                readFasta(line);
            }
        }
        else if (ARG.getStringArg("format") == "phylip") {
            ifstream fin(ARG.getStringArg("input")), fin2(ARG.getStringArg("input"));
            string line;
            while (readPhylip(fin, fin2));
        }
        else {
            cerr << "Failed to parse format named '" << ARG.getStringArg("format") << "'\n";
            exit(1);
        }
        tripInit.nSpecies = names.size();
    }

    Workflow(const Workflow &WF, const vector<int> &ind2species, size_t pos, size_t nSite, size_t offset){
        meta.initialize();
        vector<string> species2ind(WF.names.size());
        vector<int> ind2ind;
        int nInd = ind2species.size();
        for (int iInd = 0; iInd < nInd; iInd++){
            ind2ind.push_back(iInd);
            int iSpecies = ind2species[iInd];
            string name = WF.names[iSpecies] + "_" + to_string(iInd);
            names.push_back(name);
            name2id[name] = iInd;
            if (species2ind[iSpecies] == "") species2ind[iSpecies] = name;
            else species2ind[iSpecies] = string("(") + species2ind[iSpecies] + "," + name + ")";
        }
        TreeTokenizer TK(WF.guide);
        string token;
        while((token = TK()) != ""){
            if (TreeTokenizer::KEYWORDS.find(token[0]) != string::npos) guide += token;
            else guide += species2ind[WF.name2id.at(token)];
        }
        meta.nRounds = 0;
        meta.guideTree = guide;
        for (int t = 0; t < tripInit.nThreads; t++){
            size_t iSite = t * nSite / tripInit.nThreads;
            size_t new_nSite = (t + 1) * nSite / tripInit.nThreads - t * nSite / tripInit.nThreads;
            for (int iInd = 0; iInd < nInd; iInd++){
                for (int diSite = 0; diSite < new_nSite; diSite++){
                    tripInit.seq.append(WF.tripInit.seq.get(pos + iInd * offset + iSite + diSite));
                }
            }
            formatGeneNaive(ind2ind, iSite * nInd, new_nSite, new_nSite);
        }
        tripInit.nSpecies = names.size();
        tripInit.nThreads = meta.nThreads;
    }

    BinaryTree getGenetree(){
        LOG.enabled = false;
        pair<score_t, string> result = meta.run();
        LOG.enabled = true;
        return BinaryTree(result.second, name2id);
    }
};

int main(int argc, char** argv){
    ARG.setProgramName("asterisk-hky", "Accurate Species Tree EstimatoR from Independent Site Kernals");
    ARG.addStringArg('f', "format", "fasta", "Input file type, fasta: one fasta file for the whole alignment, list: a txt file containing a list of FASTA files, phylip: a phylip file for the whole alignment", true);
    ARG.addStringArg('m', "mutation", "", "Substitution rate file from Iqtree if assumming heterogeneous rates", true);
    ARG.addIntArg('d', "diskcover", 1, "The number of replicates in the disk covering method", true);
    ARG.addIntArg(0, "chunk", 10000, "The chunk size of each local region for parameter estimation");
    ARG.addIntArg(0, "iteration", 1, "The number of iterations in the iterative method");
    ARG.addFlag('I', "quick", "Set the iteration number to 1", [&]() {
		ARG.getIntArg("iteration") = 1;
	}, true);

    int iteration = 1;
    pair<score_t, string> result;
    vector<string> names;
    unordered_map<string, int> name2id;
    {
        LOG << "Iteration " << iteration << ":" << endl;
        Workflow WF(argc, argv);
        if (ARG.getStringArg("guide") != ""){
            WF.meta.nRounds = 0;
            WF.meta.nSample = 0;
        }
        if (iteration < ARG.getIntArg("iteration")) WF.meta.outputFile = "<standard output>";
        LOG << "#Base: " << WF.meta.tripInit.seq.len() << endl;
        result = WF.meta.run();
        names = WF.names;
        name2id = WF.name2id;
	}
    while (iteration < ARG.getIntArg("iteration")){
        iteration++;
        Workflow WF(argc, argv, result.second, names, name2id);
        if (iteration < ARG.getIntArg("iteration")) WF.meta.outputFile = "<standard output>";
        LOG << "Iteration " << iteration << ":" << endl;
        result = WF.meta.run();
    }
    LOG << "Score: " << (double) result.first << endl;
	return 0;
}

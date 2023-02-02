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

    void formatGene(const vector<int> &ind2species, size_t pos, size_t nSite, size_t offset) {
        int nInd = ind2species.size(), nSpecies = names.size(), nKernal = nSite, nRep = 0;
        TripartitionInitializer::Gene::Initializer gene(nInd, nSpecies, nSite, nKernal, nRep);
        buildGeneSeq(gene, ind2species, pos, nSite, offset);
        tripInit.genes.emplace_back(gene);
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

    BinaryTree getGenetree(){
        LOG.enabled = false;
        pair<score_t, string> result = meta.run();
        LOG.enabled = true;
        return BinaryTree(result.second, name2id);
    }
};

int main(int argc, char** argv){
    ARG.setProgramName("master-site", "Massive-scale Alignment-based Species Tree EstimatoR (Site)");
    ARG.addStringArg('f', "format", "fasta", "Input file type, fasta: one fasta file for the whole alignment, list: a txt file containing a list of FASTA files, phylip: a phylip file for the whole alignment", true);
    ARG.addStringArg('m', "mutation", "", "Substitution rate file from Iqtree if assumming heterogeneous rates", true);
    ARG.addIntArg('d', "diskcover", 1, "The number of replicates in the disk covering method", true);
    ARG.addIntArg(0, "chunk", 10000, "The chunk size of each local region for parameter estimation");
    
    Workflow WF(argc, argv);
    LOG << "#Base: " << WF.meta.tripInit.seq.len() << endl;
    LOG << "Score: " << (double) WF.meta.run().first << endl;
	return 0;
}

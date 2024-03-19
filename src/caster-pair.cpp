#define DRIVER_VERSION "3"

/* CHANGE LOG
 * 3: Option for no smart-pairing
 * 2: Updating input file parser
 * 1: Modified the logic in parsing FASTA names 
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
#define BLOCK_BOOTSTRAP

//#define LARGE_DATA
#ifdef LARGE_DATA
typedef long double score_t;
typedef long long count_t;
#else
typedef double score_t;
typedef int count_t;
#endif

#include "sequence-pair.hpp"
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

    void buildGeneSeq(TripartitionInitializer::Gene::Initializer &gene, const vector<int> &ind2species, size_t pos, size_t nSite, size_t offset){
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
        double total = cnt[0] + cnt[1] + cnt[2] + cnt[3];
		double p = (cnt[0] + 0.5 * (cnt[1] + cnt[2])) / total;
        gene.pq = p * (1 - p);
    }

    void formatGene(const vector<int> &ind2species, size_t pos, size_t nSite, size_t offset) {
        int nInd = ind2species.size(), nSpecies = names.size(), nKernal = nSite, nRep = 0;
        TripartitionInitializer::Gene::Initializer gene(nInd, nSpecies, nSite, nKernal, nRep);
        buildGeneSeq(gene, ind2species, pos, nSite, offset);
        tripInit.genes.emplace_back(gene);
    }

    void read(const string &file, const string &fileFormat, const string &seqFormat) {
        AlignmentParser AP(file, fileFormat, seqFormat), AP2(file, fileFormat, seqFormat);
		while (AP.nextAlignment()){
            vector<pair<long long, long long> > sitepair, sitepair2;
            AP2.nextAlignment();
            long long nSites = AP.getLength();
            {
                int d = ARG.getIntArg("pairdist");
                vector<array<unsigned short, 4> > freq;
                freq.resize(nSites);
                while (AP.nextSeq()) {
                    addName(meta.mappedname(AP.getName()));
                    string seq = AP.getSeq();
                    for (size_t i = 0; i < nSites; i++) {
                        switch (seq[i]) {
                            case 'A': freq[i][0]++; break;
                            case 'C': freq[i][1]++; break;
                            case 'G': freq[i][2]++; break;
                            case 'T': freq[i][3]++; break;
                        }
                    }
                }
                long long category[4] = {-(1 << 30), -(1 << 30), -(1 << 30), -(1 << 30)};
                for (long long i = 0; i < nSites; i++) {
                    int cnt = 0;
                    bool singleton = false;
                    for (int j = 0; j < 4; j++) {
                        if (freq[i][j] > 0) cnt++;
                        if (freq[i][j] == 1) singleton = true;
                    }
                    if (ARG.getIntArg("pairdist") > 0){
                        if (cnt == 2) { if (category[0] + d >= i) sitepair.push_back({category[0], i}); category[0] = i; }
                        if (cnt == 3) { if (category[1] + d >= i) sitepair.push_back({category[1], i}); category[1] = i; }
                        if (cnt == 4 && singleton) { if (category[2] + d >= i) sitepair.push_back({category[2], i}); category[2] = i; }
                        if (cnt == 4 && !singleton) { if (category[3] + d >= i) sitepair2.push_back({category[3], i}); category[3] = i; }
                    }
                    else {
                        if (category[0] >= 0) sitepair.push_back({category[0], i}); category[0] = i;
                    }
                }
            }
            {
                long long pos = tripInit.seq.len();
                long long offset = 3 * sitepair.size() + sitepair2.size();
                vector<int> ind2species;
                while (AP2.nextSeq()) {
                    ind2species.push_back(name2id[meta.mappedname(AP2.getName())]);
                    string seq = AP2.getSeq();
                    for (const pair<long long, long long> &e: sitepair) tripInit.seq.append<true, true, false, false>(seq[e.first], seq[e.second]);
                    for (const pair<long long, long long> &e: sitepair) tripInit.seq.append<true, false, true, false>(seq[e.first], seq[e.second]);
                    for (const pair<long long, long long> &e: sitepair) tripInit.seq.append<true, false, false, true>(seq[e.first], seq[e.second]);
                    for (const pair<long long, long long> &e: sitepair2) tripInit.seq.append<true, false, true, false>(seq[e.first], seq[e.second]);
                }
                for (int x = 0; x < 3; x++) {
                    long long len = sitepair.size();
                    long long nChunk = (len + ARG.getIntArg("chunk") - 1) / ARG.getIntArg("chunk");
                    for (int i = 0; i < nChunk; i++) {
                        size_t s = i * len / nChunk, t = (i + 1) * len / nChunk;
                        formatGene(ind2species, pos + s, t - s, offset);
                    }
                    pos += sitepair.size();
                }
                {
                    long long len = sitepair2.size();
                    long long nChunk = (len + ARG.getIntArg("chunk") - 1) / ARG.getIntArg("chunk");
                    for (int i = 0; i < nChunk; i++) {
                        long long s = i * len / nChunk, t = (i + 1) * len / nChunk;
                        formatGene(ind2species, pos + s, t - s, offset);
                    }
                }
            }
        }
    }

    Workflow(int argc, char** argv){
        meta.initialize(argc, argv);
        if (ARG.getStringArg("root") != "") addName(ARG.getStringArg("root"));
        init();
    }

    void init(){
        tripInit.nThreads = meta.nThreads;
        
        string fileFormat = ARG.getStringArg("format");
        string seqFormat = (ARG.getIntArg("datatype")) ? "AA2NA" : "NA";
        if (fileFormat != "auto" && fileFormat != "fasta" && fileFormat != "phylip" && fileFormat != "list"){
            cerr << "Failed to parse format named '" << fileFormat << "'\n";
			exit(1);
        }
		read(ARG.getStringArg("input"), fileFormat, seqFormat);
        tripInit.nSpecies = names.size();
    }
};

int main(int argc, char** argv){
    ARG.setProgramName("caster-pair", "Coalescence-aware Alignment-based Species Tree EstimatoR (Pair)");
    ARG.addStringArg('f', "format", "auto", "Input file type, fasta: one fasta file for the whole alignment, list: a txt file containing a list of FASTA files, phylip: a phylip file for the whole alignment, auto (default): detect format automatically", true);
    ARG.addStringArg('m', "mutation", "", "Substitution rate file from Iqtree if assumming heterogeneous rates", true);
    ARG.addIntArg('d', "diskcover", 1, "The number of replicates in the disk covering method", true);
    ARG.addIntArg(0, "chunk", 10000, "The chunk size of each local region for parameter estimation");
	ARG.addIntArg(0, "pairdist", 20, "The distance for pairing sites (0 for strict neighbor pairing)");
	ARG.addIntArg(0, "datatype", 0, "0 (default): nucleotides, 1: amino acids");
    
    Workflow WF(argc, argv);
    LOG << "#Base: " << WF.meta.tripInit.seq.len() << endl;
    auto res = WF.meta.run();
    LOG << "Score: " << (double) res.first << endl;
	return 0;
}

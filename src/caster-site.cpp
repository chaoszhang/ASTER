#define DRIVER_VERSION "5"

/* CHANGE LOG
 * 5: Final score normalization & Updating Eqfreq
 * 4: Updating Eqfreq calculation
 * 3: Updating input file parser
 * 2: Adding branch length functionality
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

//#define CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH
#ifndef CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH
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

    void buildGeneSeq(TripartitionInitializer::Gene::Initializer &gene, const vector<int> &ind2species, size_t pos, size_t nSite, size_t offset, const array<double, 4> &eqfreq){
        int nInd = ind2species.size();
        for (int iInd = 0; iInd < nInd; iInd++) {
            size_t pSeq = pos + iInd * offset;
            gene.species2ind[ind2species[iInd]].push_back(iInd);
            gene.ind2seq[iInd] = pSeq;
        }
        for (int j = 0; j < 4; j++) {
            gene.pi[j] = eqfreq[j];
        }
    }

    void formatGene(const vector<int> &ind2species, size_t pos, size_t nSite, size_t offset, const array<double, 4> &eqfreq) {
        int nInd = ind2species.size(), nSpecies = names.size(), nKernal = nSite, nRep = 0;
        TripartitionInitializer::Gene::Initializer gene(nInd, nSpecies, nSite, nKernal, nRep);
        buildGeneSeq(gene, ind2species, pos, nSite, offset, eqfreq);
        tripInit.genes.emplace_back(gene);
    }

    void read(const string &file, const string &fileFormat, const string &seqFormat) {
        AlignmentParser AP(file, fileFormat, seqFormat), AP2(file, fileFormat, seqFormat);
        while (AP.nextAlignment()){
            AP2.nextAlignment();
            long long nSites = AP.getLength();
            long long nChunk = (nSites + ARG.getIntArg("chunk") - 1) / ARG.getIntArg("chunk");
            vector<vector<long long> > sites(nChunk);
            vector<array<double, 4> > eqfreq;
            long long offset = 0;
            {
                vector<array<unsigned short, 4> > freq;
                freq.resize(AP.getLength());
                while (AP.nextSeq()) {
                    addName(meta.mappedname(AP.getName()));
                    string seq = AP.getSeq();
                    for (size_t i = 0; i < seq.size(); i++) {
                        switch (seq[i]) {
                            case 'A': freq[i][0]++; break;
                            case 'C': freq[i][1]++; break;
                            case 'G': freq[i][2]++; break;
                            case 'T': freq[i][3]++; break;
                        }
                    }
                }
                for (long long i = 0; i < nChunk; i++) {
                    long long s = i * nSites / nChunk, t = (i + 1) * nSites / nChunk;
                    long long sumFreq[4] = {};
                    for (long long j = s; j < t; j++) {
                        for (int k = 0; k < 4; k++) {
                            sumFreq[k] += freq[j][k];
                        }
                        #ifdef CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH
                        sites[i].push_back(j);
                        #else
                        if (freq[j][0] + freq[j][2] >= 2 && freq[j][1] + freq[j][3] >= 2) sites[i].push_back(j);
                        #endif
                    }
                    offset += sites[i].size();
                    double total = sumFreq[0] + sumFreq[1] + sumFreq[2] + sumFreq[3];
                    eqfreq.push_back({ sumFreq[0] / total, sumFreq[1] / total, sumFreq[2] / total, sumFreq[3] / total });
                }
            }
            {
                vector<int> ind2species;
                size_t pos = tripInit.seq.len();
                while (AP2.nextSeq()) {
                    ind2species.push_back(name2id[meta.mappedname(AP2.getName())]);
                    string seq = AP2.getSeq();
                    for (long long i = 0; i < nChunk; i++) {
                        for (long long j : sites[i]) tripInit.seq.append(seq[j]);
                    }
                }
                for (int i = 0; i < nChunk; i++) {
                    formatGene(ind2species, pos, sites[i].size(), offset, eqfreq[i]);
                    pos += sites[i].size();
                }
            }
        }
    }

    Workflow(int argc, char** argv){
        //string mappingFile;
        meta.initialize(argc, argv);
        if (ARG.getStringArg("root") != "") addName(ARG.getStringArg("root"));
        init();
    }

    void init(){
        tripInit.nThreads = meta.nThreads;

        string fileFormat = ARG.getStringArg("format");
        string seqFormat = (ARG.getIntArg("ambiguity")) ? "ambiguity" : "NA";
        if (fileFormat != "auto" && fileFormat != "fasta" && fileFormat != "phylip" && fileFormat != "list"){
            cerr << "Failed to parse format named '" << fileFormat << "'\n";
			exit(1);
        }
		read(ARG.getStringArg("input"), fileFormat, seqFormat);
        tripInit.nSpecies = names.size();
    }
};

int main(int argc, char** argv){
    ARG.setProgramName("caster-site", "Coalescence-aware Alignment-based Species Tree EstimatoR (Site)");
    ARG.addStringArg('f', "format", "auto", "Input file type, fasta: one fasta file for the whole alignment, list: a txt file containing a list of FASTA files, phylip: a phylip file for the whole alignment, auto (default): detect format automatically", true);
    ARG.addIntArg(0, "ambiguity", 0, "0 (default): ambiguity codes are treated as N, 1: ambiguity codes are treated as diploid unphased sites");
    ARG.addIntArg(0, "chunk", 10000, "The chunk size of each local region for parameter estimation");
    
    #ifdef CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH
    cerr << "Warning: This version can be much slower and require much more memory than the regular version. You might want to compute the correct topology first and add branch lengths with this version on fix topology.\n";
    #endif

    Workflow WF(argc, argv);
    LOG << "#Base: " << WF.meta.tripInit.seq.len() << endl;
    auto res = WF.meta.run();
    LOG << "Normalized score: " << (double) res.first / 4 << endl;
	return 0;
}

#define DRIVER_VERSION "6"

/* CHANGE LOG
 * 6: Bug fix
 * 5: Final score normalization & Switching default objective
 * 4: Option for objective functions
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
    }

    void formatGene(const vector<int> &ind2species, size_t pos, size_t nSite, size_t offset, double p) {
        int nInd = ind2species.size(), nSpecies = names.size(), nKernal = nSite, nRep = 0;
        TripartitionInitializer::Gene::Initializer gene(nInd, nSpecies, nSite, nKernal, nRep);
        buildGeneSeq(gene, ind2species, pos, nSite, offset);
        gene.pq = p * (1 - p);
        tripInit.genes.emplace_back(gene);
    }

    void read(const string &file, const string &fileFormat, const string &seqFormat) {
        AlignmentParser AP(file, fileFormat, seqFormat), AP2(file, fileFormat, seqFormat);
        while (AP.nextAlignment()){
            vector<bool> invariants;
            AP2.nextAlignment();
            long long nSites = AP.getLength();
            long long nChunk = (nSites + ARG.getIntArg("chunk") - 1) / ARG.getIntArg("chunk");
            long long d = ARG.getIntArg("pairdist");
			int objective = ARG.getIntArg("objective");
            vector<vector<pair<long long, long long> > > sitepairs(nChunk), sitepairsRY(nChunk);
            vector<double> piAC, piAG, piAT;
            vector<bool> useRY;
			int n = 0;
            long long offset = 0;
            {
                vector<array<unsigned short, 4> > freq;
                freq.resize(nSites);
                while (AP.nextSeq()) {
                    addName(meta.mappedname(AP.getName()));
                    string seq = AP.getSeq();
					n++;
                    for (long long i = 0; i < nSites; i++) {
                        switch (seq[i]) {
                            case 'A': freq[i][0]++; break;
                            case 'C': freq[i][1]++; break;
                            case 'G': freq[i][2]++; break;
                            case 'T': freq[i][3]++; break;
                        }
                    }
                }
                for (long long i = 0; i < nChunk; i++){
                    long long s = i * nSites / nChunk, t = (i + 1) * nSites / nChunk;
                    long long prev[7] = {-1 << 30, -1 << 30, -1 << 30, -1 << 30, -1 << 30, -1 << 30, -1 << 30};
                    long long sumFreq[4] = {};
					vector<int> cnts(t - s);
					vector<bool> singletons(t - s);
					long long cnt4 = 0, cnt1234 = 0;
                    for (long long j = s; j < t; j++){
                        int cnt = 0;
                        bool singleton = false;
                        for (int k = 0; k < 4; k++) {
                            if (freq[j][k] > 0) cnt++;
                            if (freq[j][k] == 1) singleton = true;
                            sumFreq[k] += freq[j][k];
                        }
                        if (cnt > 2) singleton = false;
						if (cnt >= 4) cnt4++;
						if (cnt >= 1) cnt1234++;
                        if (cnt == 4 && !singleton) cnt = 5;
                        cnts[j - s] = cnt;
						singletons[j - s] = singleton;
                        if (cnt == 0) continue;
                        if (objective == 3 && n >= 20 && d > 0){
                            long long pj = prev[cnt];
                            if (pj + d >= j && !(cnt == 1 && singletons[pj - s]) && !(cnts[pj - s] == 1 && singleton) && !(cnt == 1 && cnts[pj - s] == 1)) {
                                sitepairsRY[i].emplace_back(pj, j);
                                if (cnt < 5 && cnts[pj - s] < 5) sitepairs[i].emplace_back(pj, j);
                            }
                            // prev[cnt - 1] = prev[cnt] = prev[cnt + 1] = j;
							if (cnt > 2) prev[cnt] = j;
							else prev[1] = prev[2] = j;
                        }
                        else {
                            if (j != s && !(cnt == 1 && singletons[j - 1 - s]) && !(cnts[j - 1 - s] == 1 && singleton) && !(cnt == 1 && cnts[j - 1 - s] == 1)){
								sitepairsRY[i].emplace_back(j - 1, j);
								sitepairs[i].emplace_back(j - 1, j);
							}
                        }
                    }
					double r = (n - 1) / 60.0;
					if (objective == 1 || (objective == 3 && n < 20 && cnt4 > cnt1234 * r * r * r)) sitepairs[i].clear();
                    double totalFreq = sumFreq[0] + sumFreq[1] + sumFreq[2] + sumFreq[3];
                    piAC.push_back((sumFreq[0] + sumFreq[1]) / totalFreq);
                    piAG.push_back((sumFreq[0] + sumFreq[2]) / totalFreq);
                    piAT.push_back((sumFreq[0] + sumFreq[3]) / totalFreq);
                    offset += sitepairsRY[i].size() + 2 * sitepairs[i].size();
                }
            }

            long long pos = tripInit.seq.len();
            vector<int> ind2species;
            while (AP2.nextSeq()) {
                ind2species.push_back(name2id[meta.mappedname(AP2.getName())]);
                string seq = AP2.getSeq();
                for (long long i = 0; i < nChunk; i++){
                    for (const pair<long long, long long> &e: sitepairsRY[i]) tripInit.seq.append<true, false, true, false>(seq[e.first], seq[e.second]);
                    for (const pair<long long, long long> &e: sitepairs[i]) tripInit.seq.append<true, true, false, false>(seq[e.first], seq[e.second]);
                    for (const pair<long long, long long> &e: sitepairs[i]) tripInit.seq.append<true, false, false, true>(seq[e.first], seq[e.second]);
                }
            }
            for (int i = 0; i < nChunk; i++) {
                if (sitepairsRY[i].size()){
					formatGene(ind2species, pos, sitepairsRY[i].size(), offset, piAG[i]);
					pos += sitepairsRY[i].size();
                }
				if (sitepairs[i].size()){
                    formatGene(ind2species, pos, sitepairs[i].size(), offset, piAC[i]);
                    pos += sitepairs[i].size();
                    formatGene(ind2species, pos, sitepairs[i].size(), offset, piAT[i]);
                    pos += sitepairs[i].size();
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
        if (fileFormat != "auto" && fileFormat != "fasta" && fileFormat != "phylip" && fileFormat != "list"){
            cerr << "Failed to parse format named '" << fileFormat << "'\n";
            exit(1);
        }
        read(ARG.getStringArg("input"), fileFormat, "NA");
        tripInit.nSpecies = names.size();
    }
};

int main(int argc, char** argv){
    ARG.setProgramName("caster-pair", "Coalescence-aware Alignment-based Species Tree EstimatoR (Pair)");
    ARG.addStringArg('f', "format", "auto", "Input file type, fasta: one fasta file for the whole alignment, list: a txt file containing a list of FASTA files, phylip: a phylip file for the whole alignment, auto: detect format automatically", true);
    ARG.addIntArg(0, "chunk", 10000, "The chunk size of each local region for parameter estimation");
    ARG.addIntArg(0, "pairdist", 20, "The distance for pairing sites (0 for strict neighbor pairing)");
    ARG.addIntArg(0, "objective", 1, "Objective function, 1:RY only, 2: RY+WS+KM, 3: auto-select");

    Workflow WF(argc, argv);
    LOG << "#Base: " << WF.meta.tripInit.seq.len() << endl;
    auto res = WF.meta.run();
    LOG << "Normalized score: " << (double) res.first / 4 << endl;
    return 0;
}
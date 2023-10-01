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

string recode(const string &AA){
	string NA;
	for (const char c: AA){
		switch(c){
			case 'C': case 'c':
			case 'M': case 'm':
			case 'I': case 'i':
			case 'L': case 'l':
			case 'V': case 'v': NA.push_back('A'); break;
			case 'D': case 'd':
			case 'E': case 'e':
			case 'Q': case 'q':
			case 'N': case 'n':
			case 'H': case 'h':
			case 'R': case 'r':
			case 'K': case 'k': NA.push_back('T'); break;
			case 'S': case 's':
			case 'T': case 't':
			case 'A': case 'a':
			case 'G': case 'g':
			case 'P': case 'p': NA.push_back('C'); break;
			case 'W': case 'w':
			case 'Y': case 'y':
			case 'F': case 'f': NA.push_back('G'); break;
			default: NA.push_back('N');
		}
	}
	return NA;
}

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

    void readFasta(const string file) {
        cerr << "Processing " << file << "... \n";
        vector<pair<long long, long long> > sitepair, sitepair2;
		long long nTaxa, nSites;
        {
			int d = ARG.getIntArg("pairdist");
            vector<array<unsigned short, 4> > freq;
            ifstream fin(file);
            string name;
            fin >> name;
            while (name != "") {
                addName(meta.mappedname(name.substr(1)));
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
				if (ARG.getIntArg("datatype") == 1) seq = recode(seq);
                for (int i = 0; i < seq.size(); i++) {
                    switch (seq[i]) {
                        case 'A': case 'a': freq[i][0]++; break;
                        case 'C': case 'c': freq[i][1]++; break;
                        case 'G': case 'g': freq[i][2]++; break;
                        case 'T': case 't': freq[i][3]++; break;
						case 'U': case 'u': freq[i][3]++; break;
                    }
                }
            }
			nSites = freq.size();
			long long category[4] = {-(1 << 30), -(1 << 30), -(1 << 30), -(1 << 30)};
            for (long long i = 0; i < nSites; i++) {
                int cnt = 0;
				bool singleton = false;
                for (int j = 0; j < 4; j++) {
                    if (freq[i][j] > 0) cnt++;
					if (freq[i][j] == 1) singleton = true;
                }
                if (cnt == 2) { if (category[0] + d >= i) sitepair.push_back({category[0], i}); category[0] = i; }
				if (cnt == 3) { if (category[1] + d >= i) sitepair.push_back({category[1], i}); category[1] = i; }
				if (cnt == 4 && singleton) { if (category[2] + d >= i) sitepair.push_back({category[2], i}); category[2] = i; }
				if (cnt == 4 && !singleton) { if (category[3] + d >= i) sitepair2.push_back({category[3], i}); category[3] = i; }
            }
        }
        {
            ifstream fin(file);
            string name;
            fin >> name;
            vector<int> ind2species;
            long long pos = tripInit.seq.len();
			long long offset = 3 * sitepair.size() + sitepair2.size();
            while (name != "") {
                ind2species.push_back(name2id[meta.mappedname(name.substr(1))]);
                name = "";
                string line, seq;
                while (fin >> line) {
                    if (line[0] == '>') { name = line; break; }
                    seq += line;
                }
				if (ARG.getIntArg("datatype") == 1) seq = recode(seq);
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

    bool readPhylip(ifstream &fin, ifstream &fin2){
        long long nTaxa, nSites;
        if (!(fin >> nTaxa)) return false;
        fin >> nSites;
		vector<pair<long long, long long> > sitepair, sitepair2;
        {
			int d = ARG.getIntArg("pairdist");
            vector<array<unsigned short, 4> > freq(nSites);
            for (int i = 0; i < nTaxa; i++){
                string name, seq;
                fin >> name >> seq;
                if (seq.size() != nSites) { cerr << "The input is ill-formated."; exit(0); }
                addName(meta.mappedname(name));
                for (int j = 0; j < nSites; j++) {
                    switch (seq[j]) {
                        case 'A': case 'a': freq[j][0]++; break;
                        case 'C': case 'c': freq[j][1]++; break;
                        case 'G': case 'g': freq[j][2]++; break;
                        case 'T': case 't': freq[j][3]++; break;
						case 'U': case 'u': freq[j][3]++; break;
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
                if (cnt == 2) { if (category[0] + d >= i) sitepair.push_back({category[0], i}); category[0] = i; }
				if (cnt == 3) { if (category[1] + d >= i) sitepair.push_back({category[1], i}); category[1] = i; }
				if (cnt == 4 && singleton) { if (category[2] + d >= i) sitepair.push_back({category[2], i}); category[2] = i; }
				if (cnt == 4 && !singleton) { if (category[3] + d >= i) sitepair2.push_back({category[3], i}); category[3] = i; }
            }
        }
        {
            fin2 >> nTaxa >> nSites;
            long long pos = tripInit.seq.len();
			long long offset = 3 * sitepair.size() + sitepair2.size();
            vector<int> ind2species;
            for (int i = 0; i < nTaxa; i++){
                string name, seq;
                fin2 >> name >> seq;
                ind2species.push_back(name2id[meta.mappedname(name)]);
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
        return true;
    }

    Workflow(int argc, char** argv){
        //string mappingFile;
        meta.initialize(argc, argv);
        if (ARG.getStringArg("root") != "") addName(ARG.getStringArg("root"));
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
    ARG.setProgramName("caster-pair", "Coalescence-aware Alignment-based Species Tree EstimatoR (Pair)");
    ARG.addStringArg('f', "format", "fasta", "Input file type, fasta: one fasta file for the whole alignment, list: a txt file containing a list of FASTA files, phylip: a phylip file for the whole alignment", true);
    ARG.addStringArg('m', "mutation", "", "Substitution rate file from Iqtree if assumming heterogeneous rates", true);
    ARG.addIntArg('d', "diskcover", 1, "The number of replicates in the disk covering method", true);
    ARG.addIntArg(0, "chunk", 10000, "The chunk size of each local region for parameter estimation");
	ARG.addIntArg(0, "pairdist", 20, "The distance for pairing sites");
	ARG.addIntArg(0, "datatype", 0, "0 (default): nucleotides, 1: amino acids");
    
    Workflow WF(argc, argv);
    LOG << "#Base: " << WF.meta.tripInit.seq.len() << endl;
    auto res = WF.meta.run();
    LOG << "Score: " << (double) res.first << endl;
	return 0;
}

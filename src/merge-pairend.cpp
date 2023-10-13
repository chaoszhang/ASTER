#define DRIVER_VERSION "0"

#include<iostream>
#include<fstream>
#include<cstdio>
#include<cmath>
#include<cstdlib>
#include<cstring>
#include<string>

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

unsigned char CHAR2BITS[256] = {}, CHAR2MASK[256] = {}, RC[256] = {}, UPPER[256] = {};
unsigned char QUALITY2ASCII[95] = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";
unsigned char ASCII2QUALITY[256] = {};

const string STRING_REVERSE_COMPLEMENT(const string &seq){
    string res;
    for (int i = seq.size() - 1; i >= 0; i--){
        res += RC[seq[i]];
    }
    return res;
}

const string STRING_REVERSE(const string &seq){
    string res;
    for (int i = seq.size() - 1; i >= 0; i--){
        res += seq[i];
    }
    return res;
}

const string STRING_UPPER(const string &seq){
    string res;
    for (char c: seq){
        res += UPPER[c];
    }
    return res;
}

int MATCH_LENGTH(const string &seqL, const string &seqR, int offset, int mismatch = 1){
    int res = 0;
    for (int i = 0; i + offset < seqL.size() && i < seqR.size(); i++){
        if (seqL[i + offset] == seqR[i]) res++;
        else {
            mismatch--;
            if (mismatch == -1) return res;
        }
    }
    return res;
}

int BEST_MATCH(const string &seqL, const string &seqR, int minmatch = 20, int mismatch = 1){
    int bestmatch = minmatch - 1, bestoffset = -1;
    for (int i = 0; i < seqL.size(); i++){
        int curmatch = MATCH_LENGTH(seqL, seqR, i, mismatch);
        if (curmatch > bestmatch){
            bestmatch = curmatch;
            bestoffset = i;
        }
    }
    return bestoffset;
}

void MERGE(string &seq, string &qc, const string &seqL, const string &seqR, const string &qcL, const string &qcR, int offset){
    seq = "";
    qc = "";
    int i = 0;
    for (; i < offset; i++){
        seq += seqL[i];
        qc += qcL[i];
    }
    for (; i < seqL.size(); i++){
        if (seqL[i] == seqR[i - offset]){
            seq += seqL[i];
            qc += QUALITY2ASCII[min(93, ASCII2QUALITY[qcL[i]] + ASCII2QUALITY[qcR[i - offset]])];
        }
        else if (qcL[i] > qcR[i - offset]){
            seq += seqL[i];
            qc += QUALITY2ASCII[ASCII2QUALITY[qcL[i]] - ASCII2QUALITY[qcR[i - offset]]];
        }
        else {
            seq += seqR[i - offset];
            qc += QUALITY2ASCII[ASCII2QUALITY[qcR[i - offset]] - ASCII2QUALITY[qcL[i]]];
        }
    }
    for (; i - offset < seqR.size(); i++){
        seq += seqR[i - offset];
        qc += qcR[i - offset];
    }
}

int main(int argc, char** argv){
    if (argc != 4){
        cerr << "How to use: \nmerge-pairend FORWARD_STRAND.fastq REVERSE_STRAND.fastq MERGED_OUTPUT.fastq\n";
        exit(0);
    }

    CHAR2BITS['a'] = 0; CHAR2BITS['c'] = 1; CHAR2BITS['g'] = 2; CHAR2BITS['t'] = 3; CHAR2BITS['u'] = 3;
    CHAR2BITS['A'] = 0; CHAR2BITS['C'] = 1; CHAR2BITS['G'] = 2; CHAR2BITS['T'] = 3; CHAR2BITS['U'] = 3;

    for (size_t i = 0; i < 256; i++) CHAR2MASK[i] = 1;
    CHAR2MASK['a'] = 0; CHAR2MASK['c'] = 0; CHAR2MASK['g'] = 0; CHAR2MASK['t'] = 0; CHAR2MASK['u'] = 0;
    CHAR2MASK['A'] = 0; CHAR2MASK['C'] = 0; CHAR2MASK['G'] = 0; CHAR2MASK['T'] = 0; CHAR2MASK['U'] = 0;

    RC['a'] = 't'; RC['c'] = 'g'; RC['g'] = 'c'; RC['t'] = 'a'; RC['u'] = 'a'; 
    RC['A'] = 'T'; RC['C'] = 'G'; RC['G'] = 'C'; RC['T'] = 'A'; RC['U'] = 'A'; 

    for (size_t i = 0; i < 256; i++) UPPER[i] = 'N';
    UPPER['a'] = 'A'; UPPER['c'] = 'C'; UPPER['g'] = 'G'; UPPER['t'] = 'T'; UPPER['u'] = 'T';
    UPPER['A'] = 'A'; UPPER['C'] = 'C'; UPPER['G'] = 'G'; UPPER['T'] = 'T'; UPPER['U'] = 'T';

    for (size_t i = 0; i < 93; i++) ASCII2QUALITY[QUALITY2ASCII[i]] = i;
    
    ifstream fin1(argv[1]);
    ifstream fin2(argv[2]);
    ofstream fout(argv[3]);

    string f1line1, f1line2, f1line3, f1line4;
    string f2line1, f2line2, f2line3, f2line4;
    string seq1, seq2, qc1, qc2;

    while (getline(fin1, f1line1)){
        getline(fin1, f1line2);
        getline(fin1, f1line3);
        getline(fin1, f1line4);
        getline(fin2, f2line1);
        getline(fin2, f2line2);
        getline(fin2, f2line3);
        getline(fin2, f2line4);
        seq1 = STRING_UPPER(f1line2);
        seq2 = STRING_REVERSE_COMPLEMENT(STRING_UPPER(f2line2));
        qc1 = f1line4;
        qc2 = STRING_REVERSE(f2line4);
        int offset = BEST_MATCH(seq1, seq2);
        if (offset == -1){
            fout << f1line1 << endl;
            fout << f1line2 << endl;
            fout << f1line3 << endl;
            fout << f1line4 << endl;
            fout << f2line1 << endl;
            fout << f2line2 << endl;
            fout << f2line3 << endl;
            fout << f2line4 << endl;
        }
        else {
            string seq, qc;
            MERGE(seq, qc, seq1, seq2, qc1, qc2, offset);
            fout << "@" << endl;
            fout << seq << endl;
            fout << "+" << endl;
            fout << qc << endl;
        }
    }
    return 0;
}

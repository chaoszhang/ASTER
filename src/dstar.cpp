#define DRIVER_VERSION "1"

#include<iostream>
#include<fstream>
#include<sstream>
#include<unordered_map>
#include<cstdlib>
#include<cstring>
#include<vector>
#include<array>

using namespace std;

struct DataType16{
    typedef unsigned short FreqType;
    typedef double EqFreqType;
    typedef double ScoreType;
    typedef long long CounterType;
};

template<typename DataType> class DStarQuadrupartitionScorer{
public:
    typedef typename DataType::FreqType FreqType;
    typedef typename DataType::EqFreqType EqFreqType;
    typedef typename DataType::ScoreType ScoreType;
    typedef typename DataType::CounterType CounterType;

    struct Block{
        array<vector<FreqType>, 4> cnt0, cnt1, cnt2, cnt3;
        vector<array<EqFreqType, 4> > pi;
		int windowSize;
    };
    
private:
    inline static CounterType quadXXYY(CounterType x0, CounterType x1, CounterType x2, CounterType x3, CounterType y0, CounterType y1, CounterType y2, CounterType y3){
	    return x0 * x1 * y2 * y3 + y0 * y1 * x2 * x3;
    }

    static ScoreType scoreSite(int pos, const array<vector<FreqType>, 4> &cnt0, const array<vector<FreqType>, 4> &cnt1,
            const array<vector<FreqType>, 4> &cnt2, const array<vector<FreqType>, 4> &cnt3, const array<EqFreqType, 4> &pi){
        const EqFreqType A = pi[0], C = pi[1], G = pi[2], T = pi[3];
        const EqFreqType R = A + G, Y = C + T, R2 = A * A + G * G, Y2 = C * C + T * T;
        const FreqType a0 = cnt0[0][pos], c0 = cnt0[1][pos], g0 = cnt0[2][pos], t0 = cnt0[3][pos], r0 = a0 + g0, y0 = c0 + t0;
        const FreqType a1 = cnt1[0][pos], c1 = cnt1[1][pos], g1 = cnt1[2][pos], t1 = cnt1[3][pos], r1 = a1 + g1, y1 = c1 + t1;
        const FreqType a2 = cnt2[0][pos], c2 = cnt2[1][pos], g2 = cnt2[2][pos], t2 = cnt2[3][pos], r2 = a2 + g2, y2 = c2 + t2;
        const FreqType a3 = cnt3[0][pos], c3 = cnt3[1][pos], g3 = cnt3[2][pos], t3 = cnt3[3][pos], r3 = a3 + g3, y3 = c3 + t3;

        const CounterType rryy = quadXXYY(r0, r1, r2, r3, y0, y1, y2, y3);

        const CounterType aayy = quadXXYY(a0, a1, a2, a3, y0, y1, y2, y3);
        const CounterType ggyy = quadXXYY(g0, g1, g2, g3, y0, y1, y2, y3);
        const CounterType rrcc = quadXXYY(r0, r1, r2, r3, c0, c1, c2, c3);
        const CounterType rrtt = quadXXYY(r0, r1, r2, r3, t0, t1, t2, t3);
        
        const CounterType aacc = quadXXYY(a0, a1, a2, a3, c0, c1, c2, c3);
        const CounterType aatt = quadXXYY(a0, a1, a2, a3, t0, t1, t2, t3);
        const CounterType ggcc = quadXXYY(g0, g1, g2, g3, c0, c1, c2, c3);
        const CounterType ggtt = quadXXYY(g0, g1, g2, g3, t0, t1, t2, t3);
        
        return rryy * R2 * Y2 - (aayy + ggyy) * (R * R) * Y2 - (rrcc + rrtt) * R2 * (Y * Y)
            + (aacc + aatt + ggcc + ggtt) * (R * R) * (Y * Y);
    }

public:
    static ScoreType scoreInterval(int start, int end, const array<vector<FreqType>, 4> &cnt0, const array<vector<FreqType>, 4> &cnt1,
            const array<vector<FreqType>, 4> &cnt2, const array<vector<FreqType>, 4> &cnt3, const array<EqFreqType, 4> &pi){
        ScoreType res = 0;
        for (int i = start; i < end; i++) res += scoreSite(i, cnt0, cnt1, cnt2, cnt3, pi);
        return res;
    }
	
    static vector<ScoreType> dstar(int windowSize, const array<vector<FreqType>, 4> &cnt0, const array<vector<FreqType>, 4> &cnt1,
            const array<vector<FreqType>, 4> &cnt2, const array<vector<FreqType>, 4> &cnt3, const vector<array<EqFreqType, 4> > &pi){
        vector<ScoreType> res;
        for (int i = 0; i < cnt0[0].size(); i += windowSize){
			int j = (i + windowSize < cnt0[0].size()) ? i + windowSize : cnt0[0].size();
            res.push_back(scoreInterval(i, j, cnt0, cnt1, cnt2, cnt3, pi[i / windowSize]));
        }
        return res;
    }
	
	static CounterType quartetCnt(int start, int end, const array<vector<FreqType>, 4> &cnt0, const array<vector<FreqType>, 4> &cnt1,
            const array<vector<FreqType>, 4> &cnt2, const array<vector<FreqType>, 4> &cnt3){
		CounterType res = 0;
        for (int i = start; i < end; i++){
			CounterType s0 = cnt0[0][i] + cnt0[1][i] + cnt0[2][i] + cnt0[3][i];
			CounterType s1 = cnt1[0][i] + cnt1[1][i] + cnt1[2][i] + cnt1[3][i];
			CounterType s2 = cnt2[0][i] + cnt2[1][i] + cnt2[2][i] + cnt2[3][i];
			CounterType s3 = cnt3[0][i] + cnt3[1][i] + cnt3[2][i] + cnt3[3][i];
			res += s0 * s1 * s2 * s3;
		}
        return res;
	}
	
	static vector<CounterType> dstarQuartetCnt(int windowSize, const array<vector<FreqType>, 4> &cnt0, const array<vector<FreqType>, 4> &cnt1,
            const array<vector<FreqType>, 4> &cnt2, const array<vector<FreqType>, 4> &cnt3){
        vector<CounterType> res;
        for (int i = 0; i < cnt0[0].size(); i += windowSize){
			int j = (i + windowSize < cnt0[0].size()) ? i + windowSize : cnt0[0].size();
            res.push_back(quartetCnt(i, j, cnt0, cnt1, cnt2, cnt3));
        }
        return res;
    }

    static Block parseFreqs(const array<vector<FreqType>, 4> &f1, const array<vector<FreqType>, 4> &f2, 
            const array<vector<FreqType>, 4> &f3, const array<vector<FreqType>, 4> &f4, int start, int end, int windowSize){
        Block res;
        res.windowSize = windowSize;
        res.pi.resize((end - start + windowSize - 1) / windowSize);
        array<const array<vector<FreqType>, 4>*, 4> lst = {&f1, &f2, &f3, &f4};
        array<array<vector<FreqType>, 4>*, 4> cntlst = {&res.cnt0, &res.cnt1, &res.cnt2, &res.cnt3};
        for (int i = 0; i < 4; i++){
            const array<vector<FreqType>, 4> &f = *(lst[i]);
            array<vector<FreqType>, 4> &cnt = *(cntlst[i]);
            for (int k = 0; k < 4; k++) {
                for (int j = start; j < end; j++){
                    cnt[k].push_back(f[k][j]);
                    res.pi[(j - start) / windowSize][k] += f[k][j];
                }
            }
        }
        for (int i = 0; i < res.pi.size(); i++){
            EqFreqType sum = res.pi[i][0] + res.pi[i][1] + res.pi[i][2] + res.pi[i][3];
            for (int k = 0; k < 4; k++) res.pi[i][k] = (sum == 0) ? 0.25 : res.pi[i][k] / sum;
        }
        return res;
    }

    static string multiind(string input, string mapping = "", int intervalSize = 1000000, int windowSize = 10000, bool header = true)
    {
        string name[4];
        unordered_map<string, int> name2id;
        unordered_map<string, int> partname2id;
        if (mapping != ""){
            ifstream fmap(mapping);
            string idname, partname;
            while(fmap >> idname){
                fmap >> partname;
                if (!partname2id.count(partname)) {
                    if (partname != "-"){
                        name[partname2id.size()] = partname;
                        partname2id[partname] = partname2id.size();
                    }
                    else partname2id[partname] = -1;
                }
                name2id[idname] = partname2id[partname];
            }
        }
        ifstream fin(input);
        ostringstream fout;
        string line;
        int id, pos;
        array<array<vector<FreqType>, 4>, 4> freq;
        while (getline(fin, line)){
            if (line[0] == '>'){
                if (!name2id.count(line.substr(1))) {
                    name[partname2id.size()] = line.substr(1);
                    partname2id[line.substr(1)] = partname2id.size();
                    name2id[line.substr(1)] = partname2id[line.substr(1)];
                }
                id = name2id[line.substr(1)];
                pos = 0;
            }
            else if (id != -1){
                for (int j = 0; j < line.size(); j++){
                    for (int k = 0; k < 4; k++){
                        if (pos + j >= freq[id][k].size()) freq[id][k].push_back(0); 
                    }
                    freq[id][0][pos + j] += (line[j] == 'A' || line[j] == 'a');
                    freq[id][1][pos + j] += (line[j] == 'C' || line[j] == 'c');
                    freq[id][2][pos + j] += (line[j] == 'G' || line[j] == 'g');
                    freq[id][3][pos + j] += (line[j] == 'T' || line[j] == 't'); 
                }
                pos += line.size();
            }
        }
        if (header) fout << "file\tpos\tc*ABBA\tc*BABA\tc*AABB\tD*\tQuartetCnt\n";
        else cerr << "file\tpos\tc*ABBA\tc*BABA\tc*AABB\tD*\tQuartetCnt\n";
        
        double total1 = 0, total2 = 0, total3 = 0;
        for (int pos = 0; pos < freq[0][0].size(); pos += intervalSize){
            int end = (pos + intervalSize < freq[0][0].size()) ? pos + intervalSize : freq[0][0].size();
            Block data = parseFreqs(freq[0], freq[1], freq[2], freq[3], pos, end, windowSize);
            vector<ScoreType> topology1 = dstar(windowSize, data.cnt0, data.cnt3, data.cnt1, data.cnt2, data.pi);
            vector<ScoreType> topology2 = dstar(windowSize, data.cnt1, data.cnt3, data.cnt0, data.cnt2, data.pi);
            vector<ScoreType> topology3 = dstar(windowSize, data.cnt2, data.cnt3, data.cnt0, data.cnt1, data.pi);
            vector<CounterType> quartetCnt = dstarQuartetCnt(windowSize, data.cnt0, data.cnt1, data.cnt2, data.cnt3);
            double sum1 = 0, sum2 = 0, sum3 = 0, qcnt = 0;
            for (int i = 0; i < topology1.size(); i++){
                sum1 += topology1[i];
                sum2 += topology2[i];
                sum3 += topology3[i];
                qcnt += quartetCnt[i];
            }
            fout << input << "\t" << pos << "\t" << sum1 << "\t" << sum2 << "\t" << sum3 << "\t" << (sum1 - sum2) / (sum1 + sum2 + sum3) << "\t" << qcnt << endl;
            total1 += sum1;
            total2 += sum2;
            total3 += sum3;
        }
        
        cerr << "(((" << name[0] << "," << name[1] << ")," << name[2] << ")," << name[3] << ");\n"; 
        cerr << "c*ABBA = " << total1 << " (" << name[1] << " and " << name[2] << ")\n";
        cerr << "c*BABA = " << total2 << " (" << name[0] << " and " << name[2] << ")\n";
        cerr << "c*AABB = " << total3 << " (" << name[0] << " and " << name[1] << ")\n";
        cerr << "D* = (c*ABBA - c*BABA) / (c*ABBA + c*BABA + c*AABB) = " << (total1 - total2) / (total1 + total2 + total3) << "\n";
        
        return fout.str();
    }
};

const string HELP = R"V0G0N(D* Statistic Sliding Window Tool
dstar FASTA_FILE [ MAPPING_FILE WINDOW_SIZE ]

FASTA_FILE: input file, currently only supporting FASTA format
MAPPING_FILE: a file mapping input sequences into four clusters or - (see format below, default: -)
WINDOW_SIZE: ideally a multiple of 10000 (default: 10000)

example:
dstar input.fasta - 1000000
dstar input.fasta mapping.txt

mapping file format (exactly four clusters):
seq_name1	P1
seq_name2	P2
seq_name3	P3
seq_name4	P3
seq_name5	Po
seq_name6	Po

default mapping file format (-):
seq_name1	P1
seq_name2	P2
seq_name3	P3
seq_name4	Po
)V0G0N";

int main(int argc, char *argv[])
{
	if (argc == 1 || argv[1][0] == '-'){
		cerr << HELP;
		return 0;
	}
	
	string fasta = argv[1];
	string mapping = (argc > 2) ? argv[2] : "-";
	if (mapping == "-") mapping = "";
	int size = (argc > 3) ? stoi(argv[3]) : 10000;
	
    cout << DStarQuadrupartitionScorer<DataType16>::multiind(fasta, mapping, size);
    return 0;
}
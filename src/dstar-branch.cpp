#define DRIVER_VERSION "1"

#include<iostream>
#include<fstream>
#include<sstream>
#include<unordered_map>
#include<cstdlib>
#include<cstring>
#include<vector>
#include<array>
#include<thread>

using namespace std;

struct Tree{
    struct Node{
        bool isLeaf = true;
        Node *lc = nullptr, *rc = nullptr;
        string label;

        Node(const string &nw, int &i){
            if (nw[i] == '('){
                isLeaf = false;
                lc = new Node(nw, ++i);
                while (nw[i] != ',') i++;
                rc = new Node(nw, ++i);
                while (nw[i] != ')') i++;
            }
            else{
                int j = i;
                while (nw[i] != ':' && nw[i] != ',' && nw[i] != ')') i++;
                label = nw.substr(j, i-j); 
            }
        }

        void print(ostream& fout){
            if (!isLeaf){
                fout << '(';
                lc->print(fout);
                fout << ',';
                rc->print(fout);
                fout << ')';
            }
            fout << label;
        }

        void addLabel(int &cnt){
            if (isLeaf) return;
            lc->addLabel(cnt);
            label = string("clade") + to_string(++cnt);
            rc->addLabel(cnt);
        }

        ~Node(){
            if (isLeaf){
                delete lc;
                delete rc;
            }
        }
    } *rootNode, *ingroupNode, *outgroupNode;

    struct Instruction{
        struct Unit{
            Node* src = nullptr;
            double abbaS = 0, babaS = 0, bbaaS = 0;
            double abbaM = 0, babaM = 0, bbaaM = 0;

            Unit operator+ (const Unit &o) const{
                Unit res;
                res.abbaS = abbaS + o.abbaS;
                res.babaS = babaS + o.babaS;
                res.bbaaS = bbaaS + o.bbaaS;
                res.abbaM = abbaM + o.abbaM;
                res.babaM = babaM + o.babaM;
                res.bbaaM = bbaaM + o.bbaaM;
                return res;
            }
        };

        Node *tgt1, *tgt2, *src, *ingroup;
        bool abnormal = false, neighbor = false;
        vector<Unit> output;

        void print(ostream &fout, const string &prefix){
            for (int i = 0; i < output.size(); i++){
                if (abnormal){
                    fout << prefix << "\t" << output[i].src->label << "\t" << tgt1->label << "\t" << output[i].abbaS << "\t" << output[i].babaS << "\t" << output[i].bbaaS 
                        << "\t" << (output[i].abbaS - output[i].babaS) / (output[i].abbaS + output[i].babaS + output[i].bbaaS) << "\t" << "single" << "\t" << "undirected" << endl;
                    fout << prefix << "\t" << output[i].src->label << "\t" << tgt1->label << "\t" << output[i].abbaM << "\t" << output[i].babaM << "\t" << output[i].bbaaM 
                        << "\t" << (output[i].abbaM - output[i].babaM) / (output[i].abbaM + output[i].babaM + output[i].bbaaM) << "\t" << "multiple" << "\t" << "undirected" << endl;
                }
                else {
                    string direction = (neighbor && i == output.size() - 1) ? "directed" : "undirected";
                    fout << prefix << "\t" << output[i].src->label << "\t" << tgt1->label << "\t" << output[i].abbaS << "\t" << output[i].babaS << "\t" << output[i].bbaaS 
                        << "\t" << (output[i].abbaS - output[i].babaS) / (output[i].abbaS + output[i].babaS + output[i].bbaaS) << "\t" << "single" << "\t" << direction << endl;
                    fout << prefix << "\t" << output[i].src->label << "\t" << tgt2->label << "\t" << output[i].babaS << "\t" << output[i].abbaS << "\t" << output[i].bbaaS 
                        << "\t" << (output[i].babaS - output[i].abbaS) / (output[i].abbaS + output[i].babaS + output[i].bbaaS) << "\t" << "single" << "\t" << direction << endl;
                    fout << prefix << "\t" << output[i].src->label << "\t" << tgt1->label << "\t" << output[i].abbaM << "\t" << output[i].babaM << "\t" << output[i].bbaaM 
                        << "\t" << (output[i].abbaM - output[i].babaM) / (output[i].abbaM + output[i].babaM + output[i].bbaaM) << "\t" << "multiple" << "\t" << direction << endl;
                    fout << prefix << "\t" << output[i].src->label << "\t" << tgt2->label << "\t" << output[i].babaM << "\t" << output[i].abbaM << "\t" << output[i].bbaaM 
                        << "\t" << (output[i].babaM - output[i].abbaM) / (output[i].abbaM + output[i].babaM + output[i].bbaaM) << "\t" << "multiple" << "\t" << direction << endl;
                }
            }
        }
    };

    Tree(const string &nw){
        int i = 0;
        rootNode = new Node(nw, i);
        if (rootNode->lc->isLeaf){
            ingroupNode = rootNode->rc;
            outgroupNode = rootNode->lc;
        }
        else{
            ingroupNode = rootNode->lc;
            outgroupNode = rootNode->rc;
        }
        int cnt = 0;
        ingroupNode->lc->addLabel(cnt);
        ingroupNode->rc->addLabel(cnt);
    }

    static Instruction createInstruction(Node *ingroup, Node *tgt1, Node *tgt2, Node *src, bool neighbor, bool abnormal){
        Instruction ins;
        ins.ingroup = ingroup;
        ins.tgt1 = tgt1;
        ins.tgt2 = tgt2;
        ins.src = src;
        ins.neighbor = neighbor;
        ins.abnormal = abnormal;
        return ins;
    }

    static void createJobTgt(vector<Instruction> &jobs, Node *ingroup, Node *tgt, Node *src, bool neighbor){
        if (tgt->isLeaf) return;
        jobs.push_back(createInstruction(ingroup, tgt->lc, tgt->rc, src, neighbor, false));
        createJobTgt(jobs, ingroup, tgt->lc, src, false);
        createJobTgt(jobs, ingroup, tgt->rc, src, false);
    }

    static void createJobIngroup(vector<Instruction> &jobs, Node *ingroup){
        if (ingroup->isLeaf) return;
        createJobTgt(jobs, ingroup, ingroup->lc, ingroup->rc, true);
        createJobTgt(jobs, ingroup, ingroup->rc, ingroup->lc, true);
        if (!ingroup->lc->isLeaf){
            jobs.push_back(createInstruction(ingroup, ingroup->rc, ingroup->lc->lc, ingroup->lc->rc, true, true));
            jobs.push_back(createInstruction(ingroup, ingroup->rc, ingroup->lc->rc, ingroup->lc->lc, true, true));
        }
        if (!ingroup->rc->isLeaf){
            jobs.push_back(createInstruction(ingroup, ingroup->lc, ingroup->rc->lc, ingroup->rc->rc, true, true));
            jobs.push_back(createInstruction(ingroup, ingroup->lc, ingroup->rc->rc, ingroup->rc->lc, true, true));
        }
        createJobIngroup(jobs, ingroup->lc);
        createJobIngroup(jobs, ingroup->rc);
    }

    vector<Instruction> createJob(){
        vector<Instruction> jobs;
        createJobIngroup(jobs, ingroupNode);
        return jobs;
    }

    void print(ostream& fout){
        rootNode->print(fout);
        fout << ";\n";
    }

    ~Tree(){
        delete rootNode;
    }
};

struct DataType16{
    typedef unsigned short FreqType;
    typedef double EqFreqType;
    typedef double ScoreType;
    typedef unsigned long long CounterType;
};

struct DataType8{
    typedef unsigned char FreqType;
    typedef double EqFreqType;
    typedef double ScoreType;
    typedef unsigned int CounterType;
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

int WINDOW_SIZE = 10000;

template<typename FreqType> array<vector<FreqType>, 4>& operator+=(array<vector<FreqType>, 4>&a, const array<vector<bool>, 4> &b){
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < b[i].size(); j++){
            a[i][j] += b[i][j];
        }
    }
    return a;
}

struct Algorithm{
    typedef DataType8 DataType;
    typedef DStarQuadrupartitionScorer<DataType> Scorer;

    vector<array<vector<bool>, 4> > seqs;
    vector<string> names;
    unordered_map<string, int> name2id;
    int len, nWindow;
    array<vector<DataType::FreqType>, 4> outgroup;
    vector<array<DataType::EqFreqType, 4> > pi;

    static array<vector<bool>, 4> string2seq(const string &s){
        array<vector<bool>, 4> res;
        for (char c: s) res[0].push_back(c == 'A' || c == 'a');
        for (char c: s) res[1].push_back(c == 'C' || c == 'c');
        for (char c: s) res[2].push_back(c == 'G' || c == 'g');
        for (char c: s) res[3].push_back(c == 'T' || c == 't');
        return res;
    }

    void cntReset(array<vector<DataType::FreqType>, 4> &res){
        for (int i = 0; i < 4; i++){
            res[i].clear();
            res[i].resize(nWindow * WINDOW_SIZE, 0);
        }
    }

    void convert(array<vector<DataType::FreqType>, 4> &res, const array<vector<bool>, 4> &a){
        cntReset(res);
        res += a;
    }

    void instructionCnt(Tree::Node *node, array<vector<DataType::FreqType>, 4> &cnt, Tree::Node *excluding){
        if (node == excluding) return;
        if (node->isLeaf) cnt += seqs[name2id[node->label]];
        else{
            instructionCnt(node->lc, cnt, excluding);
            instructionCnt(node->rc, cnt, excluding);
        }
    }

    int instructionScore(Tree::Node *node, array<vector<DataType::FreqType>, 4> &cnt, Tree::Instruction &ins, Scorer::Block &block){
        if (node->isLeaf) {
            convert(cnt, seqs[name2id[node->label]]);
            Tree::Instruction::Unit unit;
            for (int i = 0; i < nWindow; i++){
                unit.abbaS += Scorer::scoreInterval(i * WINDOW_SIZE, (i + 1) * WINDOW_SIZE, block.cnt1, block.cnt2, block.cnt0, outgroup, pi[i]);
                unit.babaS += Scorer::scoreInterval(i * WINDOW_SIZE, (i + 1) * WINDOW_SIZE, block.cnt0, block.cnt2, block.cnt1, outgroup, pi[i]);
                unit.bbaaS += Scorer::scoreInterval(i * WINDOW_SIZE, (i + 1) * WINDOW_SIZE, block.cnt0, block.cnt1, block.cnt2, outgroup, pi[i]);
                unit.abbaM += Scorer::scoreInterval(i * WINDOW_SIZE, (i + 1) * WINDOW_SIZE, block.cnt1, block.cnt2, block.cnt0, block.cnt3, pi[i]);
                unit.babaM += Scorer::scoreInterval(i * WINDOW_SIZE, (i + 1) * WINDOW_SIZE, block.cnt0, block.cnt2, block.cnt1, block.cnt3, pi[i]);
                unit.bbaaM += Scorer::scoreInterval(i * WINDOW_SIZE, (i + 1) * WINDOW_SIZE, block.cnt0, block.cnt1, block.cnt2, block.cnt3, pi[i]);
            }
            ins.output.push_back(unit);
        }
        else {
            int lc = instructionScore(node->lc, cnt, ins, block);
            int rc = instructionScore(node->rc, cnt, ins, block);
            ins.output.push_back(ins.output[lc] + ins.output[rc]);
        }
        ins.output.back().src = node;
        return ins.output.size() - 1;
    }

    void runInstruction(Tree &tree, Tree::Instruction &ins){
        Scorer::Block block;
        cntReset(block.cnt0);
        cntReset(block.cnt1);
        cntReset(block.cnt2);
        cntReset(block.cnt3);
        block.pi = pi;
        block.windowSize = WINDOW_SIZE;
        instructionCnt(tree.rootNode, block.cnt3, ins.ingroup);
        instructionCnt(ins.tgt2, block.cnt0, nullptr);
        if (ins.abnormal) {
            instructionCnt(ins.tgt1, block.cnt2, nullptr);
            instructionScore(ins.src, block.cnt1, ins, block);
        }
        else {
            instructionCnt(ins.tgt1, block.cnt1, nullptr);
            instructionScore(ins.src, block.cnt2, ins, block);
        }
    }

    void threadJobs(Tree &tree, vector<Tree::Instruction> &jobs, int L, int R){
        for (int i = L; i < R; i++) runInstruction(tree, jobs[i]);
    }

    Algorithm(const string &faFile, const string &nwFile, const string &tsvFile, int nThread){
        ifstream ffa(faFile), fnw(nwFile);
        string nw, line, tempName, tempSeq;
        getline(fnw, nw);
        Tree tree(nw);
        tree.print(cout);
        while(getline(ffa, line)){
            if (line[0] == '>'){
                if (tempSeq.size()){
                    name2id[tempName] = names.size();
                    names.push_back(tempName);
                    seqs.push_back(string2seq(tempSeq));
                }
                tempName = line.substr(1);
                tempSeq.clear();
            }
            else tempSeq += line;
        }
        name2id[tempName] = names.size();
        names.push_back(tempName);
        seqs.push_back(string2seq(tempSeq));
        len = tempSeq.size();
        nWindow = (len - 1) / WINDOW_SIZE + 1;
        
        array<vector<DataType::FreqType>, 4> cnt;
        for (int i = 0; i < 4; i++) cnt[i].resize(nWindow * WINDOW_SIZE, 0);
        for (int i = 0; i < seqs.size(); i++) cnt += seqs[i];
        vector<array<DataType::EqFreqType, 4> > sum(nWindow);
        pi.resize(nWindow);
        for (int i = 0; i < 4; i++){
            for (int j = 0; j < len; j++){
                sum[j / WINDOW_SIZE][i] += cnt[i][j];
            }
        }
        for (int i = 0; i < nWindow; i++){
            DataType::EqFreqType total = sum[i][0] + sum[i][1] + sum[i][2] + sum[i][3];
            for (int j = 0; j < 4; j++){
                pi[i][j] = (total > 0) ? sum[i][j] / total : 0.25;
            }
        }
        convert(outgroup, seqs[name2id[tree.outgroupNode->label]]);
        
        vector<Tree::Instruction> jobs = tree.createJob();
        int nJob = jobs.size();
        vector<thread> threads;
        for (int i = 1; i < nThread; i++){
            threads.emplace_back(&Algorithm::threadJobs, this, ref(tree), ref(jobs), i * nJob / nThread, (i + 1) * nJob / nThread);
        }
        threadJobs(tree, jobs, 0, nJob / nThread);
        for (thread &t: threads) t.join();
        ofstream fout(tsvFile);
        fout << "data" << "\t" << "src" << "\t" << "tgt" << "\t" << "abba" << "\t" << "baba" << "\t" << "bbaa" << "\t" << "dstar" << "\t" << "outgroup" << "\t" << "direction" << endl;
        for (Tree::Instruction ins: jobs) ins.print(fout, faFile);
    }
};

const string HELP = R"V0G0N(D* branch
dstar-branch FASTA_FILE TREE_FILE OUTPUT_FILE [ NUM_THREADS ]

FASTA_FILE: alignment file, currently only supporting FASTA format
TREE_FILE: species tree file, currently only supporting NEWICK format and already rooted at one outgroup species
OUTPUT_TILE: output file in TSV format
NUM_THREADS: number of threads (default: 1)

example:
dstar-branch alignment.fasta speices.nw output.tsv
dstar-branch alignment.fasta speices.nw output.tsv 128
)V0G0N";

int main(int argc, char *argv[])
{
    if (argc <= 3){
		cerr << HELP;
		return 0;
	}
	
	int nThread = (argc > 4) ? stoi(argv[4]) : 1;
	
    Algorithm alg(argv[1], argv[2], argv[3], nThread);
    return 0;
}
#define DRIVER_VERSION "0"

#include<iostream>
#include<fstream>
#include<sstream>
#include<unordered_map>
#include<cstdlib>
#include<cstring>
#include<vector>
#include<array>

using namespace std;

struct Row{
	double abba, baba, bbaa;
};

struct Birow{
	double dstarbranch;
	bool directed;
};

const string HELP = R"V0G0N(Merge D* Branch Outputs & Take Diagnal Minimum
dstar-branch-merge-diag TSV_FILE1 [ TSV_FILE2 TSV_FILE3 ... ]

TSV_FILE1: output file of dstar-branch
TSV_FILE2: another output file for the same species tree

outputs:
dstar-branch.merged.tsv: merged results
dstar-branch.merged.diag.tsv: merged results taking diagnal mininum

example:
dstar-branch alignment.fasta speices.nw output.tsv
dstar-branch alignment.fasta speices.nw output.tsv 128
)V0G0N";

int main(int argc, char *argv[])
{
    if (argc == 1){
		cerr << HELP;
		return 0;
	}
	
	string data, src, tgt, abba, baba, bbaa, dstar, outgroup, direction;
	vector<Row> rows;
	unordered_map<string, unordered_map<string, Birow> > single, multiple;
	
	{
		ifstream fin(argv[1]);
		fin >> data >> src >> tgt >> abba >> baba >> bbaa >> dstar >> outgroup >> direction;
		while (fin >> data){
			Row r;
			fin >> src >> tgt >> r.abba >> r.baba >> r.bbaa >> dstar >> outgroup >> direction;
			rows.push_back(r);
		}
	}
	
	for (int i = 2; i < argc; i++){
		ifstream fin(argv[i]);
		fin >> data >> src >> tgt >> abba >> baba >> bbaa >> dstar >> outgroup >> direction;
		for (int j = 0; j < rows.size(); j++){
			Row r;
			fin >> data >> src >> tgt >> r.abba >> r.baba >> r.bbaa >> dstar >> outgroup >> direction;
			rows[j].abba += r.abba;
			rows[j].baba += r.baba;
			rows[j].bbaa += r.bbaa;
		}
	}
	
	{
		ifstream fin(argv[1]);
		ofstream fout("dstar-branch.merged.tsv");
		fin >> data >> src >> tgt >> abba >> baba >> bbaa >> dstar >> outgroup >> direction;
		fout << "src" << "\t" << "tgt" << "\t" << "abba" << "\t" << "baba" << "\t" << "bbaa" << "\t" << "dstar" << "\t" << "dstarbranch" << "\t" << "outgroup" << "\t" << "direction" << endl;
		for (int j = 0; j < rows.size(); j++){
			fin >> data >> src >> tgt >> abba >> baba >> bbaa >> dstar >> outgroup >> direction;
			Birow br;
			double rdstar = (rows[j].abba - rows[j].baba) / (rows[j].abba + rows[j].baba + rows[j].bbaa);
			br.dstarbranch = max(rdstar, 0.0);
			br.directed = (direction == "directed");
			fout << src << "\t" << tgt << "\t" << rows[j].abba << "\t" << rows[j].baba << "\t" << rows[j].bbaa << "\t" << rdstar << "\t" << br.dstarbranch << "\t" << outgroup << "\t" << direction << endl;
			if (outgroup == "single") single[src][tgt] = br;
			else multiple[src][tgt] = br;
		}
	}
	
	{
		ifstream fin(argv[1]);
		ofstream fout("dstar-branch.merged.diag.tsv");
		fin >> data >> src >> tgt >> abba >> baba >> bbaa >> dstar >> outgroup >> direction;
		fout << "src" << "\t" << "tgt" << "\t" << "dstarbranch" << "\t" << "outgroup" << "\t" << "direction" << endl;
		for (int j = 0; j < rows.size(); j++){
			fin >> data >> src >> tgt >> abba >> baba >> bbaa >> dstar >> outgroup >> direction;
			Birow br1 = (outgroup == "single") ? single[src][tgt] : multiple[src][tgt];
			Birow br2 = (outgroup == "single") ? single[tgt][src] : multiple[tgt][src];
			direction = (br1.directed && br2.directed) ? "bidirected" : "unknown";
			fout << src << "\t" << tgt << "\t" << min(br1.dstarbranch, br2.dstarbranch) << "\t" << outgroup << "\t" << direction << endl;
		}
	}
	
    return 0;
}
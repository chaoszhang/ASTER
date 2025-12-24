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
dstar-branch-merge-diag FILE_LIST

TSV_FILE1: output file of dstar-branch
TSV_FILE2: another output file for the same species tree
FILE_LIST: a list of dstar-branch tsv files in a list file, one line per each tsv file

outputs:
dstar-branch.merged.tsv: merged results
dstar-branch.merged.diag.tsv: merged results taking diagnal mininum

example:
dstar-branch-merge-diag chr1.tsv chr2.tsv chr3.tsv
dstar-branch-merge-diag chr_tsv_list.txt
)V0G0N";

int main(int argc, char *argv[])
{
    if (argc == 1){
		cerr << HELP;
		return 0;
	}
	
	bool isTsv = true;
	vector<string> files;
	string data, src, tgt, abba, baba, bbaa, dstar, outgroup, direction;
	vector<Row> rows;
	unordered_map<string, unordered_map<string, Birow> > single, multiple;
	
	if (argc == 2) {
		string line;
		ifstream fin(argv[1]);
		getline(fin, line);
		stringstream sline(line);
		if (isTsv && !(sline >> data && data == "data")) isTsv = false;
		if (isTsv && !(sline >> src && src == "src")) isTsv = false;
		if (isTsv && !(sline >> tgt && tgt == "tgt")) isTsv = false;
		if (isTsv && !(sline >> abba && abba == "abba")) isTsv = false;
		if (isTsv && !(sline >> baba && baba == "baba")) isTsv = false;
		if (isTsv && !(sline >> bbaa && bbaa == "bbaa")) isTsv = false;
		if (isTsv && !(sline >> dstar && dstar == "dstar")) isTsv = false;
		if (isTsv && !(sline >> outgroup && outgroup == "outgroup")) isTsv = false;
		if (isTsv && !(sline >> direction && direction == "direction")) isTsv = false;
	}

	if (isTsv) {
		for (int i = 1; i < argc; i++) files.push_back(argv[i]);
	}
	else {
		string line;
		ifstream fin(argv[1]);
		while(getline(fin, line)) files.push_back(line);
	}

	{
		ifstream fin(files[0]);
		fin >> data >> src >> tgt >> abba >> baba >> bbaa >> dstar >> outgroup >> direction;
		while (fin >> data){
			Row r;
			fin >> src >> tgt >> r.abba >> r.baba >> r.bbaa >> dstar >> outgroup >> direction;
			rows.push_back(r);
		}
	}
	
	for (int i = 1; i < files.size(); i++){
		ifstream fin(files[i]);
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
		ifstream fin(files[0]);
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
		ifstream fin(files[0]);
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
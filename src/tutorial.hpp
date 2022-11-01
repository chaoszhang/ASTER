#ifndef TUTORIAL
#define TUTORIAL

#include<regex>
#include<fstream>

class MDGenerator{
const string APRO_UNIQUE_INTRO = R"V0G0N(# Accurate Species Tree ALgorithm for PaRalogs and Orthologs (ASTRAL-Pro2)
ASTRAL-Pro stands for ASTRAL for PaRalogs and Orthologs. ASTRAL is a tool for estimating an unrooted species tree given a set of unrooted gene trees and is statistically consistent under the multi-species coalescent model (and thus is useful for handling incomplete lineage sorting, i.e., ILS). ASTRAL-pro extends ASTRAL to allow multi-copy genes. ASTRAL-pro finds the species tree that has the maximum number of shared induced quartet tree equivalent classes with the set of gene trees, subject to the constraint that the set of bipartitions in the species tree comes from a predefined set of bipartitions. Please see the paper below for the definition of the PL-quartet scores, which is what ASTRAL-Pro optimizes. We refer to the tool both as A-Pro and ASTRAL-Pro. 

`ASTRAL-Pro2` re-implements [ASTRAL-Pro](https://github.com/chaoszhang/A-pro) in an equally accurate yet **faster**, and **easier to install** and **lower memory consumption** way.

## Publication

[1] Chao Zhang, Siavash Mirarab, ASTRAL-Pro 2: ultrafast species tree reconstruction from multi-copy gene family trees, Bioinformatics, 2022, btac620, https://doi.org/10.1093/bioinformatics/btac620

[2] Chao Zhang, Celine Scornavacca, Erin K Molloy, Siavash Mirarab, ASTRAL-Pro: Quartet-Based Species-Tree Inference despite Paralogy, Molecular Biology and Evolution, Volume 37, Issue 11, November 2020, Pages 3292–3307, https://doi.org/10.1093/molbev/msaa139

### Example of usage

We obtained the species tree from muti-copy gene family trees using ASTRAL-Pro2 VERSION [1] by optimizing the objective function of ASTRAL-Pro [2].

)V0G0N";

const string ASTRAL_UNIQUE_INTRO = R"V0G0N(# Accurate Species Tree ALgorithm (wASTRAL-unweighted)
ASTRAL is a tool for estimating an unrooted species tree given a set of unrooted gene trees. ASTRAL is statistically consistent under the multi-species coalescent model (and thus is useful for handling incomplete lineage sorting, i.e., ILS). ASTRAL finds the species tree that has the maximum number of shared induced quartet trees with the set of gene trees, subject to the constraint that the set of tripartitions in the species tree comes from a predefined set of tripartitions.

`wASTRAL-unweighted` re-implements [ASTRAL](https://github.com/smirarab/ASTRAL) as a scalable alternative to ASTRAL on datasets for which ASTRAL is not suitable (e.g. large datasets, multi-individual, and gene trees with missing taxa).

As a scalable alternative to ASTRAL-III, wASTRAL-unweighted lacks of some features of ASTRAL-III (e.g. bootstrapping). You can work around by first computing optimal tree with wASTRAL-unweighted and use the wASTRAL-unweighted output tree as `-q` option to ASTRAL-III. 

## Publication

[1] Chao Zhang, Siavash Mirarab, Weighting by Gene Tree Uncertainty Improves Accuracy of Quartet-based Species Trees, Molecular Biology and Evolution, 2022, msac215, https://doi.org/10.1093/molbev/msac215

[2] Chao Zhang, Maryam Rabiee, Erfan Sayyari, and Siavash Mirarab. 2018. “ASTRAL-III: Polynomial Time Species Tree Reconstruction from Partially Resolved Gene Trees.” BMC Bioinformatics 19 (S6): 153. [doi:10.1186/s12859-018-2129-y](https://doi.org/10.1186/s12859-018-2129-y).

### Example of usage

We obtained the species tree from gene trees using wASTRAL-unweighted VERSION [1] by optimizing the objective function of ASTRAL [2].

)V0G0N";

const string SHARED_INTRO = R"V0G0N(
## Bug Reports

Contact ``aster-users@googlegroups.com`` or post on [ASTER issues page](https://github.com/chaoszhang/ASTER/issues).

# Documentations
- The rest of this TUTORIAL file
- [README/astral-pro.md](README/astral-pro.md) for ASTRAL and ASTRAL-Pro; [README/wastral.md](README/wastral.md) for weighted ASTRAL series; [README/asterisk.md](README/asterisk.md) for ASTERISK series
- Forums:
  - [User group discussions](https://groups.google.com/forum/#!forum/aster-users)
  - [ASTER issues page](https://github.com/chaoszhang/ASTER/issues)

# INSTALLATION
For most users, installing ASTER is ***very*** easy!
Download using one of two approaches:
  - You simply need to download the zip file for [Windows](https://github.com/chaoszhang/ASTER/archive/refs/heads/Windows.zip)/[MacOS](https://github.com/chaoszhang/ASTER/archive/refs/heads/MacOS.zip)/[Linux](https://github.com/chaoszhang/ASTER/archive/refs/heads/Linux.zip) and extract the contents to a folder of your choice.
  - Alternatively, you can clone the [github repository](https://github.com/chaoszhang/ASTER.git) and checkout the branch named Windows/MacOS/Linux.

Binary files should be in the `exe` folder for Windows or `bin` folder otherwise. If you are lucky, these may just work as is and you may not need to build at all.

## For Linux/Unix/WSL users
1. In terminal, `cd` into the downloaded directory and run `make`.
  - If you see `*** Installation complete! ***` then you are done!
  - If you see `Command 'g++' not found` then before rerunning `make`,
    - Debian (Ubuntu) users try
      ```
      sudo apt update
      sudo apt install g++
      ```
    - CentOS (RedHat) users try
      ```
      sudo yum update
      sudo yum install gcc-c++
      ```
    - Unix (MacOS) users should be prompted for installing `g++` and please click "install". If no prompt, try `g++`
2. Binary files should be in the `bin` folder.

## For Windows users
- Executables for x86-64 are available in `exe` folder and it is **very likely** that they already work.
- [Windows Subsystem for Linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/install) is HIGHLY recommanded if you need to install on your own! Please follow instructions in "For Linux/Unix/WSL users" section.
- To compile windows excutables:
  1. Download [MinGW](https://sourceforge.net/projects/mingw-w64/) and install ***posix*** version for your architecture (eg. x86-64)
  2. Add path to `bin` folder of MinGW to [system environment variable `PATH`](https://www.google.com/search?q=Edit+the+system+environment+variables+windows)
  3. Double click `make.bat` inside the downloaded directory

### GUI for Windows users (NEW)

Please check out our software with GUI. Simply download the [zip file](https://github.com/chaoszhang/ASTER/archive/refs/heads/Windows.zip), extract the contents, enter `exe` folder, and click `aster-gui.exe`. 
)V0G0N";

const string APRO_IO = R"V0G0N(
# INPUT
* The input gene trees are in the Newick format
* The input trees can have missing taxa, polytomies (unresolved branches), and multi-copy genes.
* When multiple genes from the same species are available, you can ask ASTRAL to force them to be together in the species tree. You can do this in two ways.
  1. You can give multiple genes from the same species the same name in the input gene trees (e.g., `((species_name_A,species_name_B),(species_name_A,species_name_C));`).
  2. OR, a mapping file needs to be provided using the `-a` option. This mapping file should have one line per genes, and each line needs to be in the following formats (e.g., for gene trees like `((gene_A1,gene_B1),(gene_A2,gene_C1));`):

```
gene_A1 species_name_A
gene_A2 species_name_A
gene_B1 species_name_B
gene_B2 species_name_B
gene_B3 species_name_B
...
```

# OUTPUT
The output in is Newick format and gives:

* the species tree topology
* branch lengths in coalescent units (only for internal branches)
* branch supports measured as [local posterior probabilities](http://mbe.oxfordjournals.org/content/early/2016/05/12/molbev.msw079.short?rss=1)
* It can also annotate branches with other quantities, such as quartet supports and localPPs for all three topologies.

The ASTRAL-Pro tree leaves the branch length of terminal branches empty. Some tools for visualization and tree editing do not like this (e.g., ape). In FigTree, if you open the tree several times, it eventually opens up (at least on our machines). In ape, if you ask it to ignore branch lengths all together, it works. In general, if your tool does not like the lack of terminal branches, you can add a dummy branch length, [as in this script](https://github.com/smirarab/global/blob/master/src/mirphyl/utils/add-bl.py).
)V0G0N";

const string ASTRAL_IO = R"V0G0N(
# INPUT
* The input trees can have missing taxa, polytomies (unresolved branches), and multiple individuals per species.
* When individuals genes from the same species are available, you can ask ASTRAL to force them to be together in the species tree. You can do this in two ways.
  1. You can give multiple individuals from the same species the same name in the input gene trees (e.g., `((species_name_A,species_name_B),(species_name_A,species_name_C));`).
  2. OR, a mapping file needs to be provided using the `-a` option. This mapping file should have one line per genes, and each line needs to be in the following formats (e.g., for gene trees like `((individual_A1,individual_B1),(individual_A2,individual_C1));`):
```
individual_A1 species_name_A
individual_A2 species_name_A
individual_B1 species_name_B
individual_B2 species_name_B
individual_B3 species_name_B
...
```

# OUTPUT
The output in is Newick format and gives:

* the species tree topology
* branch lengths in coalescent units (only for internal branches)
* branch supports measured as [local posterior probabilities](http://mbe.oxfordjournals.org/content/early/2016/05/12/molbev.msw079.short?rss=1)
* It can also annotate branches with other quantities, such as quartet supports and localPPs for all three topologies.

The ASTRAL tree leaves the branch length of terminal branches empty. Some tools for visualization and tree editing do not like this (e.g., ape). In FigTree, if you open the tree several times, it eventually opens up (at least on our machines). In ape, if you ask it to ignore branch lengths all together, it works. In general, if your tool does not like the lack of terminal branches, you can add a dummy branch length, [as in this script](https://github.com/smirarab/global/blob/master/src/mirphyl/utils/add-bl.py).
)V0G0N";

const string SHARED_EXE = R"V0G0N(
# EXECUTION
ASTER currently has no GUI. You need to run it through the command-line. In a terminal/PowerShell, go to the directory (location) where you have downloaded ASTER and issue the following command:

```
bin/PROGRAM_NAME
```

This will give you a list of options available. If you are using Windows, please replace `bin/PROGRAM_NAME` with `.\exe\PROGRAM_NAME.exe`.

To find the species tree with input from in a file called `INPUT_FILE`, use:

```
bin/PROGRAM_NAME INPUT_FILE
```
or
```
bin/PROGRAM_NAME -i INPUT_FILE
```

In the first case, INPUT_FILE is ***hard-coded*** to be the ***last argument*** for backward compatibility. 

For example if you want to run `PROGRAM_NAME` with input `example/EXAMPLE_INPUT`, then run

```
bin/PROGRAM_NAME example/EXAMPLE_INPUT
```
or
```
bin/PROGRAM_NAME -i example/EXAMPLE_INPUT
```

The results will be outputted to the standard output. To save the results in a file use the `-o OUTPUT_FILE` option before `INPUT_FILE`(**Strongly recommended**):

```
bin/PROGRAM_NAME -o OUTPUT_FILE INPUT_FILE
```
or
```
bin/PROGRAM_NAME -i INPUT_FILE -o OUTPUT_FILE
```

With `-i INPUT_FILE` option, the order does not matter anymore. For brevity, from here on we will not demonstrate `-i INPUT_FILE` cases.

To save the logs (**also recommended**), run:

```
bin/PROGRAM_NAME -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```

For example, you can run

```
bin/PROGRAM_NAME -o example/EXAMPLE_INPUT.stree example/EXAMPLE_INPUT 2>example/EXAMPLE_INPUT.log
```

ASTER supports multi-threading. To run program with 4 threads, add `-t 4` before `INPUT_FILE`:

```
bin/PROGRAM_NAME -t 4 -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```

ASTER has very good parrallel efficiency up to 64 cores when input data is large. In fact, it often experiences super-linear speedup with 16 cores or more. So feel free to use as many cores as you want.
)V0G0N";

const string APRO_UNIQUE_EXE = R"V0G0N(
By default, ASTRAL-Pro assumes multiple genes from the same species in the same input gene trees having the same name. Alternatively, a mapping file needs to be provided using the `-a` option (see INPUT section). For example,

```
bin/astral-pro -a example/multitree_genename.map example/multitree_genename.nw
```
)V0G0N";

const string ASTRAL_UNIQUE_EXE = R"V0G0N(
By default, ASTRAL assumes multiple individuals from the same species in the same input gene trees having the same name. Alternatively, a mapping file needs to be provided using the `-a` option (see INPUT section). For example,

```
bin/astral -a example/genetree.map example/genetree.nw
```

When your dataset has no more than 50 species and no more than 500 genes, you may want to run with more rounds using `-R` (see below). 
)V0G0N";

const string SHARED_ADV = R"V0G0N(
## Advanced Options

ASTER algorithm first performs `R` (4 by default) rounds of search and then repeatedly performs `S` (4 by default) rounds of subsampling and exploration until no improvement found.

```
bin/PROGRAM_NAME -r R -s S -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```

If you want to run with more rounds of placement for ensured optimality, then you can run with
```
bin/PROGRAM_NAME -r 16 -s 16 -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```
or simply
```
bin/PROGRAM_NAME -R -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```

If you want to place taxa on an existing ***fully resolved*** species tree, you can use `-c SPECIES_TREE_IN_NEWICK_FORMAT` before `INPUT_FILE`:

```
bin/PROGRAM_NAME -o OUTPUT_FILE -c SPECIES_TREE_IN_NEWICK_FORMAT INPUT_FILE
```

Specifically, you can score and annotate a ***fully resolved*** species tree containing all taxa with `-c SPECIES_TREE_IN_NEWICK_FORMAT`. If want to score a species tree or you want to place only ***one*** taxon onto the tree, you can use

```
bin/PROGRAM_NAME -r 1 -s 0 -o OUTPUT_FILE -c SPECIES_TREE_IN_NEWICK_FORMAT INPUT_FILE
```
or simply,
```
bin/PROGRAM_NAME -C -o OUTPUT_FILE -c SPECIES_TREE_IN_NEWICK_FORMAT INPUT_FILE
```

If you want to give hints by providing candidate species trees or trees similar to the species tree, you can use `-g SPECIES_TREES_IN_NEWICK_FORMAT` before `INPUT_FILE`:

```
bin/PROGRAM_NAME -o OUTPUT_FILE -g SPECIES_TREES_IN_NEWICK_FORMAT INPUT_FILE
```
)V0G0N";

const string APRO_UNIQUE_ADV = R"V0G0N(
Species tree with more than **5000** taxa may cause **overflow**. Use the following command instead:

```
bin/astral-pro_int128 -o OUTPUT_FILE INPUT_FILE
```

Add `-u 0` before `INPUT_FILE` if you want to compute species tree topology only; Add `-u 2` before `INPUT_FILE` if you support and local-PP for all three resolutions of each branch.

```
bin/astral-pro -u 0 -o OUTPUT_FILE INPUT_FILE
bin/astral-pro -u 2 -o OUTPUT_FILE INPUT_FILE
```

If you do not want to compute optimal species tree but instead just want to root and tag gene trees, you can use the following command:

```
bin/astral-pro -T -o OUTPUT_FILE INPUT_FILE
```
)V0G0N";

const string ASTRAL_UNIQUE_ADV = R"V0G0N(
Species tree with more than **5000** taxa may cause **overflow**. Use the following command instead:

```
bin/astral_int128 -o OUTPUT_FILE INPUT_FILE
```
)V0G0N";

    const string APRO = "astral-pro";
    const string ASTRAL = "astral";

    string replace(string txt, string from, string to){
        return regex_replace(txt, regex(from), to);
    }

    string uniqueIntro(){
        string rawtext = "# Accurate Species Tree EstimatoR (ASTER❋)";
        if (programName == APRO) rawtext = APRO_UNIQUE_INTRO;
        if (programName == ASTRAL) rawtext = ASTRAL_UNIQUE_INTRO;
        return replace(rawtext, "VERSION", version);
    }

    string sharedIntro(){
        return SHARED_INTRO;
    }

    string inputOutput(){
        if (programName == APRO) return APRO_IO;
        if (programName == ASTRAL) return ASTRAL_IO;

        return "";
    }

    string sharedExe(){
        if (programName == APRO) return replace(replace(SHARED_EXE, "PROGRAM_NAME", "astral-pro"), "EXAMPLE_INPUT", "multitree.nw");
        if (programName == ASTRAL) return replace(replace(SHARED_EXE, "PROGRAM_NAME", "astral"), "EXAMPLE_INPUT", "genetree.nw");

        return SHARED_EXE;
    }

    string uniqueExe(){
        if (programName == APRO) return APRO_UNIQUE_EXE;
        if (programName == ASTRAL) return ASTRAL_UNIQUE_EXE;

        return "";
    }

    string sharedAdv(){
        if (programName == APRO) return replace(replace(SHARED_ADV, "PROGRAM_NAME", "astral-pro"), "EXAMPLE_INPUT", "multitree.nw");
        if (programName == ASTRAL) return replace(replace(SHARED_ADV, "PROGRAM_NAME", "astral"), "EXAMPLE_INPUT", "genetree.nw");

        return SHARED_ADV;
    }

    string uniqueAdv(){
        if (programName == APRO) return APRO_UNIQUE_ADV;
        if (programName == ASTRAL) return ASTRAL_UNIQUE_ADV;

        return "";
    }

    string programName;

public:
    
    static string version;

    MDGenerator(string programName): programName(programName){
        string ui = uniqueIntro();
        string si = sharedIntro();
        string io = inputOutput();
        string se = sharedExe();
        string ue = uniqueExe();
        string sa = sharedAdv();
        string ua = uniqueAdv();
        ofstream fout(string("tutorial/") + programName + ".md");
        fout << (ui + si + io + se + ue + sa + ua);
    }

};

string MDGenerator::version;

#endif